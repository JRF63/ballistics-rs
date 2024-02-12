pub type FloatType = f64;
pub type Vector3 = nalgebra::base::Vector3<FloatType>;

pub type State = nalgebra::base::Vector6<FloatType>;

pub trait StateTrait {
    fn from_pos_and_vel(x0: Vector3, v0: Vector3) -> Self;

    fn pos(&self) -> Vector3;

    fn vel(&self) -> Vector3;
}

impl StateTrait for State {
    fn from_pos_and_vel(x0: Vector3, v0: Vector3) -> Self {
        State::new(x0.x, x0.y, x0.z, v0.x, v0.y, v0.z)
    }

    fn pos(&self) -> Vector3 {
        Vector3::new(self.x, self.y, self.z)
    }

    fn vel(&self) -> Vector3 {
        Vector3::new(self.w, self.a, self.b)
    }
}
