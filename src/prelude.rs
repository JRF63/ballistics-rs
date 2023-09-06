pub type FloatType = f64;
pub type Vector3 = glam::f64::DVec3;
pub const GRAVITY_ACCEL: Vector3 = Vector3::new(0.0, 0.0, -9.80665);
pub const PI: FloatType = std::f64::consts::PI;

pub trait AngleConversion {
    fn to_arcminute(self) -> Self;
}

impl AngleConversion for FloatType {
    fn to_arcminute(self) -> Self {
        self * 10800.0 / PI
    }
}
