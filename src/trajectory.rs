use crate::{
    data::lerp,
    environment::{calc_air_density, calc_speed_sound},
    solver::OdeSolver,
    state::{FloatType, State, Vec3, GRAVITY_ACCEL},
};
use std::cmp::Ordering;

pub fn calc_trajectory<F, G>(
    x0: Vec3,
    v0: Vec3,
    cd_func: F,
    wind: Vec3,
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
        let dv = vw * (-(air_density * cd_func(mach_num)) * speed) + GRAVITY_ACCEL;

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
    x0: Vec3,
    muzzle_speed: FloatType,
    cd_func: F,
    zero_range: FloatType,
    zero_elevation: FloatType,
    wind: Vec3,
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
            * Vec3::new(
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
                        drop = lerp(
                            state_old.pos.x,
                            state_old.pos.z,
                            state.pos.x,
                            state.pos.z,
                            zero_range,
                        );
                        windage = lerp(
                            state_old.pos.x,
                            state_old.pos.y,
                            state.pos.x,
                            state.pos.y,
                            zero_range,
                        );
                        true
                    }
                },
                None => false,
            }
        };

        calc_trajectory(
            x0,
            v_guess,
            cd_func,
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
