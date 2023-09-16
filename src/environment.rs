use crate::prelude::*;

/// Calculates the saturation vapor pressure of water using the Arden Buck equation.
///
/// https://en.wikipedia.org/wiki/Arden_Buck_equation
pub fn calc_saturation_vapor_pressure_water(temp: FloatType) -> FloatType {
    let temp_celsius = temp - 273.15;

    // Coefficients adjusted to give the pressure in Pascals
    let (a, b, c, d) = if temp_celsius > 0.0 {
        (611.21, 18.678, 234.5, 257.14)
    } else {
        (611.15, 23.036, 333.7, 279.82)
    };

    a * FloatType::exp((b - temp_celsius / c) * (temp_celsius / (d + temp_celsius)))
}

/// Calculates for the density of air, assuming air is an ideal gas.
pub fn calc_air_density(temp: FloatType, pressure: FloatType, rh: FloatType) -> FloatType {
    // Precomputed because Rust does not support const eval for floats.
    // Equal to:
    // let R: FloatType = 8.31446261815324; // J/(K*mol)
    // let M: FloatType = 0.0289645; // kg/mol
    // let rm: FloatType = R / M;
    let rm: FloatType = 287.0570048905812;
    let density = pressure / (rm * temp);

    // Equation 8.24 of Modern Exterior Ballistics
    // Presumably corrects for the lower molar mass of humid air
    let p_wv = calc_saturation_vapor_pressure_water(temp);
    let correction_factor = 1.0 - 0.00378 * rh * p_wv / 101325.0;

    density * correction_factor
}

/// Calculates the speed of sound, assuming air is an ideal gas.
pub fn calc_speed_sound(temp: FloatType, _pressure: FloatType, rh: FloatType) -> FloatType {
    // Precomputed because Rust does not support const eval for floats.
    // Equal to:
    // let gamma: FloatType = 1.4; // unitless
    // let R: FloatType = 8.31446261815324; // J/(K*mol)
    // let M: FloatType = 0.0289645; // kg/mol
    // let c: FloatType = (gamma * R / M).sqrt();
    let c: FloatType = 20.046940086876443;
    let speed = c * temp.sqrt();

    // Equation 8.26 of Modern Exterior Ballistics
    // Presumably corrects for the lower molar mass of humid air
    let p_wv = calc_saturation_vapor_pressure_water(temp);
    let correction_factor = 1.0 + 0.0014 * rh * p_wv / 101325.0;

    speed * correction_factor
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::almost_equal;

    const EPSILON: FloatType = 1e-3;

    #[test]
    fn test_calc_air_density() {
        let temp = 273.15;
        let pressure = 101325.0;
        let rh = 0.0;

        assert!(almost_equal(
            calc_air_density(temp, pressure, rh),
            1.2922521350727452,
            EPSILON
        ))
    }

    #[test]
    fn test_calc_speed_sound() {
        let temp = 273.15;
        let pressure = 101325.0;
        let rh = 0.0;

        assert!(almost_equal(
            calc_speed_sound(temp, pressure, rh),
            331.3207950615342,
            EPSILON
        ))
    }

    #[test]
    fn test_calc_water_vapor_pressure() {
        // https://en.wikipedia.org/wiki/Vapour_pressure_of_water
        // Units in (degC, kPa)
        let reference: &[(FloatType, FloatType)] = &[
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

        for (temp_celsius, pressure_ref) in reference {
            let pressure = calc_saturation_vapor_pressure_water(temp_celsius + 273.15);
            assert!(almost_equal(pressure, pressure_ref * 1e3, EPSILON));
        }
    }
}
