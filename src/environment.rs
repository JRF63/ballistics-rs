use uom::si::{
    f64::*,
    mass_density::kilogram_per_cubic_meter,
    pressure::pascal,
    thermodynamic_temperature::{degree_celsius, kelvin},
    velocity::meter_per_second,
};

/// Calculates the saturation vapor pressure of water using the Arden Buck equation.
///
/// https://en.wikipedia.org/wiki/Arden_Buck_equation
pub fn calc_saturation_vapor_pressure_water(temp: ThermodynamicTemperature) -> Pressure {
    let temp_celsius = temp.get::<degree_celsius>();

    // Coefficients adjusted to give the pressure in Pascals
    let (a, b, c, d) = if temp_celsius > 0.0 {
        (611.21, 18.678, 234.5, 257.14)
    } else {
        (611.15, 23.036, 333.7, 279.82)
    };

    Pressure::new::<pascal>(
        a * f64::exp((b - temp_celsius / c) * (temp_celsius / (d + temp_celsius))),
    )
}

/// Calculates for the density of air, assuming air is an ideal gas.
pub fn calc_air_density(
    temp: ThermodynamicTemperature,
    pressure: Pressure,
    rh: f64,
) -> MassDensity {
    // Precomputed because Rust does not support const eval for floats.
    // Equal to:
    // let R: f64 = 8.31446261815324; // J/(K*mol)
    // let M: f64 = 0.0289645; // kg/mol
    // let rm: f64 = R / M;
    let rm: f64 = 287.0570048905812;
    let density = pressure.get::<pascal>() / (rm * temp.get::<kelvin>());

    // Equation 8.24 of Modern Exterior Ballistics
    // Presumably corrects for the lower molar mass of humid air
    let p_wv = calc_saturation_vapor_pressure_water(temp);
    let correction_factor = 1.0 - 0.00378 * rh * p_wv.get::<pascal>() / 101325.0;

    MassDensity::new::<kilogram_per_cubic_meter>(density * correction_factor)
}

/// Calculates the speed of sound, assuming air is an ideal gas.
pub fn calc_speed_sound(temp: ThermodynamicTemperature, _pressure: Pressure, rh: f64) -> Velocity {
    // Precomputed because Rust does not support const eval for floats.
    // Equal to:
    // let gamma: f64 = 1.4; // unitless
    // let R: f64 = 8.31446261815324; // J/(K*mol)
    // let M: f64 = 0.0289645; // kg/mol
    // let c: f64 = (gamma * R / M).sqrt();
    let c: f64 = 20.046940086876443;
    let speed = c * temp.get::<kelvin>().sqrt();

    // Equation 8.26 of Modern Exterior Ballistics
    // Presumably corrects for the lower molar mass of humid air
    let p_wv = calc_saturation_vapor_pressure_water(temp);
    let correction_factor = 1.0 + 0.0014 * rh * p_wv.get::<pascal>() / 101325.0;

    Velocity::new::<meter_per_second>(speed * correction_factor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::pressure::kilopascal;

    const EPSILON: f64 = 1e-3;

    fn almost_equal(val: f64, reference: f64, epsilon: f64) -> bool {
        (val - reference).abs() < epsilon * reference.abs()
    }

    #[test]
    fn test_calc_air_density() {
        let temp = ThermodynamicTemperature::new::<degree_celsius>(0.0);
        let pressure = Pressure::new::<pascal>(101325.0);
        let rh = 0.0;

        assert!(almost_equal(
            calc_air_density(temp, pressure, rh).get::<kilogram_per_cubic_meter>(),
            1.2922521350727452,
            EPSILON
        ))
    }

    #[test]
    fn test_calc_speed_sound() {
        let temp = ThermodynamicTemperature::new::<degree_celsius>(0.0);
        let pressure = Pressure::new::<pascal>(101325.0);
        let rh = 0.0;

        assert!(almost_equal(
            calc_speed_sound(temp, pressure, rh).get::<meter_per_second>(),
            331.3207950615342,
            EPSILON
        ))
    }

    #[test]
    fn test_calc_water_vapor_pressure() {
        // https://en.wikipedia.org/wiki/Vapour_pressure_of_water
        // Units in (degC, kPa)
        let reference: &[(f64, f64)] = &[
            (0.0, 0.6113),
            (5.0, 0.8726),
            (10.0, 1.2281),
            (15.0, 1.7056),
            (20.0, 2.3388),
            (25.0, 3.1690),
            (30.0, 4.2455),
            (35.0, 5.6267),
            (40.0, 7.3814),
            (45.0, 9.5898),
            (50.0, 12.3440),
            (55.0, 15.7520),
            (60.0, 19.9320),
            (65.0, 25.0220),
            (70.0, 31.1760),
            (75.0, 38.5630),
            (80.0, 47.3730),
            (85.0, 57.8150),
            (90.0, 70.1170),
            (95.0, 84.5290),
            (100.0, 101.3200),
        ];

        for (temp, pressure_ref) in reference {
            let pressure = calc_saturation_vapor_pressure_water(ThermodynamicTemperature::new::<
                degree_celsius,
            >(*temp));
            assert!(almost_equal(
                pressure.get::<kilopascal>(),
                *pressure_ref,
                EPSILON
            ));
        }
    }
}
