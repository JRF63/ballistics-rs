use uom::si::{f64::*, pressure::hectopascal, thermodynamic_temperature::degree_celsius};

/// Calculates the vapor pressure of water using the Arden Buck equation.
pub fn calc_water_vapor_pressure(temp: ThermodynamicTemperature) -> Pressure {
    let temp_celsius = temp.get::<degree_celsius>();
    let (a, b, c, d) = if temp_celsius > 0.0 {
        (611.21, 18.678, 234.5, 257.14)
    } else {
        (611.5, 23.036, 333.7, 279.82)
    };

    Pressure::new::<hectopascal>(
        a * f64::exp((b - temp_celsius / c) * (temp_celsius / (d + temp_celsius))),
    )
}
