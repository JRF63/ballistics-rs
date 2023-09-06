use crate::{
    data::lerp,
    environment::{calc_air_density, calc_speed_sound},
    prelude::*,
    solver::OdeSolver,
    state::State,
};
use std::cmp::Ordering;

pub fn calc_trajectory<F, G>(
    x0: Vector3,
    v0: Vector3,
    drag_func: F,
    wind: Vector3,
    temp: FloatType,
    pressure: FloatType,
    rh: FloatType,
    _elevation: FloatType,
    dt: FloatType,
    t_max: FloatType,
    mut stop_eval: G,
) where
    F: Fn(FloatType) -> FloatType,
    G: FnMut(&State, FloatType) -> bool,
{
    let derivative = |state: State| -> State {
        let air_density = calc_air_density(temp, pressure, rh);
        let speed_sound = calc_speed_sound(temp, pressure, rh);

        let vw = state.vel - wind;
        let speed = vw.length();
        let mach_num = speed / speed_sound;

        let dx = state.vel;
        let dv = vw * (-(air_density * drag_func(mach_num)) * speed) + GRAVITY_ACCEL;

        State::new(dx, dv)
    };

    let y0 = State::new(x0, v0);
    let mut solver = OdeSolver::new(y0, dt);
    loop {
        let (state, time) = solver.step(derivative);
        if time > t_max || stop_eval(state, time) {
            break;
        }
    }
}

pub fn solve_initial_velocity<F>(
    x0: Vector3,
    muzzle_speed: FloatType,
    drag_func: F,
    zero_range: FloatType,
    zero_elevation: FloatType,
    wind: Vector3,
    temp: FloatType,
    pressure: FloatType,
    rh: FloatType,
    elevation: FloatType,
    dt: FloatType,
    t_max: FloatType,
) -> Option<(FloatType, FloatType)>
where
    F: Fn(FloatType) -> FloatType + Copy,
{
    const MAX_CONVERGENCE_STEPS: u32 = 100;
    const CONVERGENCE_EPSILON: FloatType = 1e-5;

    let mut ver_angle = (zero_elevation / zero_range).atan();
    let mut ver_angle_low = ver_angle - FloatType::to_radians(60.0);
    let mut ver_angle_high = ver_angle + FloatType::to_radians(60.0);
    let mut ver_angle_old = FloatType::NAN;

    let mut hor_angle = 0.0;
    let mut hor_angle_left = -FloatType::to_radians(60.0);
    let mut hor_angle_right = FloatType::to_radians(60.0);
    let mut hor_angle_old = FloatType::NAN;

    let mut ver_converged = false;
    let mut hor_converged = false;

    for _ in 0..MAX_CONVERGENCE_STEPS {
        if ver_converged && hor_converged {
            break;
        }

        if ver_angle == ver_angle_old && hor_angle == hor_angle_old {
            break;
        }

        ver_angle_old = ver_angle;
        hor_angle_old = hor_angle;

        let v_guess = muzzle_speed
            * Vector3::new(
                ver_angle.cos() * hor_angle.cos(),
                hor_angle.sin(),
                ver_angle.sin() * hor_angle.cos(),
            );

        let mut state_old = State::new(x0, v_guess);
        let mut drop = FloatType::NAN;
        let mut windage = FloatType::NAN;

        let range_reached = |state: &State, _: FloatType| -> bool {
            match state.pos.x.partial_cmp(&zero_range) {
                Some(ord) => match ord {
                    Ordering::Less => {
                        state_old = state.clone();
                        false
                    }
                    Ordering::Equal => {
                        drop = state.pos.z;
                        windage = state.pos.y;
                        true
                    }
                    Ordering::Greater => {
                        let alpha =
                            (zero_range - state_old.pos.x) / (state.pos.x - state_old.pos.x);
                        drop = lerp(state_old.pos.z, state.pos.z, alpha);
                        windage = lerp(state_old.pos.y, state.pos.y, alpha);
                        true
                    }
                },
                None => false,
            }
        };

        calc_trajectory(
            x0,
            v_guess,
            drag_func,
            wind,
            temp,
            pressure,
            rh,
            elevation,
            dt,
            t_max,
            range_reached,
        );

        // Second zero should be attained at the specified distance
        if (drop - zero_elevation).abs() < CONVERGENCE_EPSILON {
            ver_converged = true;
        } else {
            // Lost convergence, retry again with larger bounds
            if ver_converged {
                ver_converged = false;
                ver_angle_high += ver_angle_high - ver_angle;
                ver_angle_low += ver_angle_low - ver_angle;
            }

            if drop > zero_elevation {
                // Aiming too high
                ver_angle_high = ver_angle;
            } else {
                // Aiming too low
                ver_angle_low = ver_angle;
            }
        }

        if windage.abs() < CONVERGENCE_EPSILON {
            hor_converged = true;
        } else {
            if hor_converged {
                hor_converged = false;
                hor_angle_right += hor_angle_right - hor_angle;
                hor_angle_left += hor_angle_left - hor_angle;
            }

            if windage > 0.0 {
                // Aiming too far to the right
                hor_angle_right = hor_angle;
            } else {
                // Aiming too far to the left
                hor_angle_left = hor_angle;
            }
        }

        ver_angle = (ver_angle_low + ver_angle_high) / 2.0;
        hor_angle = (hor_angle_left + hor_angle_right) / 2.0;
    }

    if !(ver_converged && hor_converged) {
        None
    } else {
        Some((ver_angle, hor_angle))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{almost_equal, create_standard_drag_function, StandardDragFunction};

    const EPSILON: FloatType = 5e-2;

    const KG_PER_LB: FloatType = 0.45359237;
    const M_PER_IN: FloatType = 0.0254;
    const PASCALS_PER_MILLIBAR: FloatType = 100.0;
    const KELVIN_OFFSET: FloatType = 273.15;

    // Example 8.1 of Modern Exterior Ballistics
    #[test]
    fn test_vacuum_trajectory() {
        let muzzle_speed = 80.0;
        let firing_angle = FloatType::to_radians(2005.03 / 60.0);
        let x0 = Vector3::ZERO;
        let v0 = muzzle_speed * Vector3::new(firing_angle.cos(), 0.0, firing_angle.sin());
        let wind = Vector3::ZERO;
        let temp = 273.15;
        let pressure = 101325.0;
        let rh = 0.0;
        let elevation = 0.0;
        let dt = 0.01;
        let t_max = 10.0;

        let drag_func = |_: FloatType| -> FloatType { 0.0 };

        let stop_eval = |state: &State, t: FloatType| -> bool {
            let pos = x0 + v0 * t + GRAVITY_ACCEL / 2.0 * t * t;
            let vel = v0 + GRAVITY_ACCEL * t;

            assert!((state.pos - pos).length() < EPSILON);
            assert!((state.vel - vel).length() < EPSILON);
            false
        };

        calc_trajectory(
            x0, v0, drag_func, wind, temp, pressure, rh, elevation, dt, t_max, stop_eval,
        );
    }

    fn calc_trajectory_on_range_intervals<F>(
        x0: Vector3,
        v0: Vector3,
        drag_func: F,
        wind: Vector3,
        temp: FloatType,
        pressure: FloatType,
        rh: FloatType,
        elevation: FloatType,
        dt: FloatType,
        t_max: FloatType,
        ranges: Vec<FloatType>,
    ) -> Vec<(State, FloatType)>
    where
        F: Fn(FloatType) -> FloatType,
    {
        let mut output = Vec::with_capacity(ranges.len());

        let mut ranges_iter = ranges.iter();
        let mut range_current = *ranges_iter.next().unwrap();

        let mut state_old = State::new(x0, v0);
        let mut t_old = 0.0;

        let stop_eval = |state: &State, t: FloatType| -> bool {
            let res = loop {
                match state.pos.x.partial_cmp(&range_current) {
                    Some(ord) => match ord {
                        Ordering::Less => {
                            break false;
                        }
                        Ordering::Equal | Ordering::Greater => {
                            let mut save = |state: &State, t: FloatType| {
                                output.push((state.clone(), t));
                            };

                            if ord == Ordering::Equal {
                                save(state, t);
                            } else {
                                let alpha = (range_current - state_old.pos.x)
                                    / (state.pos.x - state_old.pos.x);
                                let state_lerp = lerp(state_old, state.clone(), alpha);
                                let t_lerp = lerp(t_old, t, alpha);
                                save(&state_lerp, t_lerp);
                            };

                            match ranges_iter.next() {
                                Some(r) => range_current = *r,
                                None => break true,
                            }
                        }
                    },
                    None => break true,
                }
            };
            state_old = state.clone();
            t_old = t;
            res
        };

        calc_trajectory(
            x0, v0, drag_func, wind, temp, pressure, rh, elevation, dt, t_max, stop_eval,
        );

        output
    }

    #[test]
    fn test_flat_fire_trajectory_no_wind() {
        // Validation data using JBM Ballistics calc
        // https://www.jbmballistics.com/cgi-bin/jbmtraj-5.1.cgi
        let reference: &[(FloatType, FloatType, FloatType, FloatType, FloatType)] = &[
            (0.0, -1.5, 0.0, 1000.0, 0.000),
            (100.0, -3.5, -3.1, 933.3, 0.104),
            (200.0, -10.0, -4.4, 869.9, 0.215),
            (300.0, -21.6, -6.3, 809.3, 0.334),
            (400.0, -39.0, -8.5, 751.2, 0.462),
            (500.0, -63.4, -11.1, 695.6, 0.600),
            (600.0, -95.7, -13.9, 642.3, 0.750),
            (700.0, -137.5, -17.1, 591.5, 0.912),
            (800.0, -190.2, -20.8, 543.3, 1.089),
            (900.0, -256.2, -24.9, 498.1, 1.281),
            (1000.0, -337.7, -29.5, 456.4, 1.491),
            (1100.0, -437.8, -34.8, 418.7, 1.720),
            (1200.0, -560.0, -40.7, 385.7, 1.969),
            (1300.0, -708.2, -47.6, 358.1, 2.238),
            (1400.0, -886.5, -55.3, 336.1, 2.527),
            (1500.0, -1099.1, -64.0, 318.8, 2.834),
            (1600.0, -1349.8, -73.7, 304.9, 3.155),
            (1700.0, -1642.2, -84.4, 293.1, 3.491),
            (1800.0, -1979.8, -96.0, 282.9, 3.839),
            (1900.0, -2366.1, -108.7, 273.7, 4.201),
            (2000.0, -2804.5, -122.4, 265.4, 4.574),
        ];

        let muzzle_speed = 1000.0;
        let bc = 0.5 * KG_PER_LB / (M_PER_IN * M_PER_IN);
        let sight_height = 1.5 * M_PER_IN;
        let wind = Vector3::ZERO;
        let temp = 25.0 + KELVIN_OFFSET; // 25 degC
        let pressure = 1013.25 * PASCALS_PER_MILLIBAR; // 1013.25 millibars
        let rh = 0.0;
        let elevation = 0.0;
        let dt = 0.01;
        let t_max = 10.0;

        let x0 = Vector3::new(0.0, 0.0, elevation - sight_height);
        let v0 = Vector3::new(muzzle_speed, 0.0, 0.0);

        let drag_func = create_standard_drag_function(StandardDragFunction::G1, bc);

        let ranges: Vec<_> = reference.iter().map(|x| x.0).collect();

        let trajectory = calc_trajectory_on_range_intervals(
            x0, v0, drag_func, wind, temp, pressure, rh, elevation, dt, t_max, ranges,
        );

        assert_eq!(reference.len(), trajectory.len());

        for ((state, t), ref_data) in trajectory.iter().zip(reference) {
            let x = state.pos.x.round();
            let z = (state.pos.z / M_PER_IN * 10.0).round() / 10.0;
            let speed = (state.vel.length() * 10.0).round() / 10.0;
            let t = (t * 1e3).round() / 1e3;

            let drop_moa = if state.pos.x == 0.0 {
                0.0
            } else {
                ((state.pos.z / state.pos.x).atan().to_arcminute() * 10.0).round() / 10.0
            };

            assert!(almost_equal(x, ref_data.0, EPSILON));
            assert!(almost_equal(z, ref_data.1, EPSILON));
            assert!(almost_equal(drop_moa, ref_data.2, EPSILON));
            assert!(almost_equal(speed, ref_data.3, EPSILON));
            assert!(almost_equal(t, ref_data.4, EPSILON));
        }
    }

    #[test]
    fn test_flat_fire_trajectory_strong_wind() {
        // Validation data using JBM Ballistics calc
        // https://www.jbmballistics.com/cgi-bin/jbmtraj-5.1.cgi
        let reference: &[(
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
        )] = &[
            (0.0, -1.5, 0.0, 0.0, 0.0, 1000.0, 0.000),
            (100.0, -3.5, -3.1, 0.7, 0.6, 933.3, 0.104),
            (200.0, -10.0, -4.4, 2.9, 1.2, 869.9, 0.215),
            (300.0, -21.6, -6.3, 6.6, 1.9, 809.3, 0.334),
            (400.0, -39.0, -8.5, 12.2, 2.7, 751.2, 0.462),
            (500.0, -63.4, -11.1, 19.7, 3.4, 695.6, 0.600),
            (600.0, -95.7, -13.9, 29.5, 4.3, 642.3, 0.750),
            (700.0, -137.5, -17.1, 41.8, 5.2, 591.5, 0.912),
            (800.0, -190.2, -20.8, 56.8, 6.2, 543.3, 1.089),
            (900.0, -256.2, -24.9, 75.0, 7.3, 498.1, 1.281),
            (1000.0, -337.7, -29.5, 96.6, 8.4, 456.4, 1.491),
            (1100.0, -437.8, -34.8, 122.0, 9.7, 418.7, 1.720),
            (1200.0, -560.0, -40.7, 151.4, 11.0, 385.7, 1.969),
            (1300.0, -708.2, -47.6, 184.7, 12.4, 358.1, 2.238),
            (1400.0, -886.5, -55.3, 221.9, 13.8, 336.1, 2.527),
            (1500.0, -1099.1, -64.0, 262.5, 15.3, 318.8, 2.834),
            (1600.0, -1349.8, -73.7, 306.1, 16.7, 304.9, 3.155),
            (1700.0, -1642.2, -84.4, 352.5, 18.1, 293.1, 3.491),
            (1800.0, -1979.9, -96.0, 401.5, 19.5, 282.9, 3.839),
            (1900.0, -2366.2, -108.7, 452.9, 20.8, 273.7, 4.201),
            (2000.0, -2804.6, -122.4, 506.7, 22.1, 265.4, 4.574),
        ];

        let muzzle_speed = 1000.0;
        let bc = 0.5 * KG_PER_LB / (M_PER_IN * M_PER_IN);
        let sight_height = 1.5 * M_PER_IN;
        let wind_speed = 5.0; // 5 m/s to the right
        let wind = Vector3::new(0.0, wind_speed, 0.0);
        let temp = 25.0 + KELVIN_OFFSET; // 25 degC
        let pressure = 1013.25 * PASCALS_PER_MILLIBAR; // 1013.25 millibars
        let rh = 0.0;
        let elevation = 0.0;
        let dt = 0.01;
        let t_max = 10.0;

        let x0 = Vector3::new(0.0, 0.0, elevation - sight_height);
        let v0 = Vector3::new(muzzle_speed, 0.0, 0.0);

        let drag_func = create_standard_drag_function(StandardDragFunction::G1, bc);

        let ranges: Vec<_> = reference.iter().map(|x| x.0).collect();

        let trajectory = calc_trajectory_on_range_intervals(
            x0, v0, drag_func, wind, temp, pressure, rh, elevation, dt, t_max, ranges,
        );

        assert_eq!(reference.len(), trajectory.len());

        for ((state, t), ref_data) in trajectory.iter().zip(reference) {
            let x = state.pos.x.round();
            let z = (state.pos.z / M_PER_IN * 10.0).round() / 10.0;
            let y = (state.pos.y / M_PER_IN * 10.0).round() / 10.0;
            let speed = (state.vel.length() * 10.0).round() / 10.0;
            let t = (t * 1e3).round() / 1e3;

            let drop_moa = if state.pos.x == 0.0 {
                0.0
            } else {
                ((state.pos.z / state.pos.x).atan().to_arcminute() * 10.0).round() / 10.0
            };

            let windage_moa = if state.pos.x == 0.0 {
                0.0
            } else {
                ((state.pos.y / state.pos.x).atan().to_arcminute() * 10.0).round() / 10.0
            };

            assert!(almost_equal(x, ref_data.0, EPSILON));
            assert!(almost_equal(z, ref_data.1, EPSILON));
            assert!(almost_equal(drop_moa, ref_data.2, EPSILON));
            assert!(almost_equal(y, ref_data.3, EPSILON));
            assert!(almost_equal(windage_moa, ref_data.4, EPSILON));
            assert!(almost_equal(speed, ref_data.5, EPSILON));
            assert!(almost_equal(t, ref_data.6, EPSILON));
        }
    }

    #[test]
    fn test_firing_angle_estimation_from_zeroing() {
        let reference: &[(
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
            FloatType,
        )] = &[
            (0.0, -1.5, 0.0, 0.0, 0.0, 1000.0, 0.000),
            (100.0, 1.5, 1.3, -0.7, -0.6, 933.3, 0.104),
            (200.0, -0.0, -0.0, 0.0, 0.0, 869.9, 0.215),
            (300.0, -6.6, -1.9, 2.4, 0.7, 809.2, 0.334),
            (400.0, -19.1, -4.2, 6.5, 1.4, 751.2, 0.462),
            (500.0, -38.4, -6.7, 12.6, 2.2, 695.6, 0.600),
            (600.0, -65.8, -9.6, 20.9, 3.0, 642.3, 0.750),
            (700.0, -102.5, -12.8, 31.8, 4.0, 591.5, 0.912),
            (800.0, -150.3, -16.4, 45.4, 5.0, 543.3, 1.089),
            (900.0, -211.2, -20.5, 62.1, 6.0, 498.1, 1.281),
            (1000.0, -287.8, -25.1, 82.3, 7.2, 456.4, 1.491),
            (1100.0, -382.9, -30.4, 106.3, 8.4, 418.7, 1.720),
            (1200.0, -500.1, -36.4, 134.2, 9.8, 385.7, 1.969),
            (1300.0, -643.3, -43.2, 166.2, 11.2, 358.1, 2.238),
            (1400.0, -816.6, -50.9, 201.9, 12.6, 336.1, 2.527),
            (1500.0, -1024.2, -59.6, 241.1, 14.0, 318.8, 2.834),
            (1600.0, -1269.9, -69.3, 283.3, 15.5, 304.9, 3.155),
            (1700.0, -1557.3, -80.0, 328.2, 16.9, 293.1, 3.491),
            (1800.0, -1890.0, -91.7, 375.7, 18.2, 282.9, 3.839),
            (1900.0, -2271.3, -104.4, 425.7, 19.6, 273.7, 4.201),
            (2000.0, -2704.6, -118.1, 478.1, 20.9, 265.4, 4.574),
        ];

        let muzzle_speed = 1000.0;
        let bc = 0.5 * KG_PER_LB / (M_PER_IN * M_PER_IN);
        let sight_height = 1.5 * M_PER_IN;
        let wind_speed = 5.0; // 5 m/s to the right
        let wind = Vector3::new(0.0, wind_speed, 0.0);
        let temp = 25.0 + KELVIN_OFFSET; // 25 degC
        let pressure = 1013.25 * PASCALS_PER_MILLIBAR; // 1013.25 millibars
        let rh = 0.0;
        let elevation = 0.0;
        let dt = 0.01;
        let t_max = 10.0;

        let zero_range = 200.0;
        let zero_elevation = 0.0;

        let x0 = Vector3::new(0.0, 0.0, elevation - sight_height);

        let drag_func = create_standard_drag_function(StandardDragFunction::G1, bc);

        let (ver_angle, hor_angle) = solve_initial_velocity(
            x0,
            muzzle_speed,
            drag_func,
            zero_range,
            zero_elevation,
            wind,
            temp,
            pressure,
            rh,
            elevation,
            dt,
            t_max,
        )
        .unwrap();

        let v0 = muzzle_speed
            * Vector3::new(
                ver_angle.cos() * hor_angle.cos(),
                hor_angle.sin(),
                ver_angle.sin() * hor_angle.cos(),
            );

        let drag_func = create_standard_drag_function(StandardDragFunction::G1, bc);

        let ranges: Vec<_> = reference.iter().map(|x| x.0).collect();

        let trajectory = calc_trajectory_on_range_intervals(
            x0, v0, drag_func, wind, temp, pressure, rh, elevation, dt, t_max, ranges,
        );

        assert_eq!(reference.len(), trajectory.len());

        for ((state, t), ref_data) in trajectory.iter().zip(reference) {
            let x = state.pos.x.round();
            let z = (state.pos.z / M_PER_IN * 10.0).round() / 10.0;
            let y = (state.pos.y / M_PER_IN * 10.0).round() / 10.0;
            let speed = (state.vel.length() * 10.0).round() / 10.0;
            let t = (t * 1e3).round() / 1e3;

            let drop_moa = if state.pos.x == 0.0 {
                0.0
            } else {
                ((state.pos.z / state.pos.x).atan().to_arcminute() * 10.0).round() / 10.0
            };

            let windage_moa = if state.pos.x == 0.0 {
                0.0
            } else {
                ((state.pos.y / state.pos.x).atan().to_arcminute() * 10.0).round() / 10.0
            };

            assert!(almost_equal(x, ref_data.0, EPSILON));
            assert!(almost_equal(z, ref_data.1, EPSILON));
            assert!(almost_equal(drop_moa, ref_data.2, EPSILON));
            assert!(almost_equal(y, ref_data.3, EPSILON));
            assert!(almost_equal(windage_moa, ref_data.4, EPSILON));
            assert!(almost_equal(speed, ref_data.5, EPSILON));
            assert!(almost_equal(t, ref_data.6, EPSILON));
        }
    }
}
