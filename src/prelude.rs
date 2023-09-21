pub use crate::{
    data::{create_standard_drag_function, cubic_hermite_interpolation, StandardDragFunction},
    solver::OdeSolver,
    state::State,
    trajectory::{calc_trajectory, solve_initial_velocity},
};

pub type FloatType = f64;
pub type Vector3 = glam::f64::DVec3;
pub(crate) const GRAVITY_ACCEL: Vector3 = Vector3::new(0.0, 0.0, -9.80665);
pub(crate) const PI: FloatType = core::f64::consts::PI;

pub trait AngleConversion {
    fn to_arcminute(self) -> Self;
}

impl AngleConversion for FloatType {
    fn to_arcminute(self) -> Self {
        self * 10800.0 / PI
    }
}
