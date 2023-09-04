use std::ops::{Add, AddAssign, Mul, MulAssign};

pub type FloatType = f64;
pub type Vec3 = glam::f64::DVec3;
pub const GRAVITY_ACCEL: Vec3 = Vec3::new(0.0, 0.0, -9.80665);

#[derive(Clone, Copy)]
pub struct State {
    pub pos: Vec3,
    pub vel: Vec3,
}

impl State {
    pub fn new(x0: Vec3, v0: Vec3) -> Self {
        Self { pos: x0, vel: v0 }
    }
}

impl Add<State> for State {
    type Output = State;

    fn add(self, rhs: State) -> Self::Output {
        State {
            pos: self.pos + rhs.pos,
            vel: self.vel + rhs.vel,
        }
    }
}

impl AddAssign<State> for State {
    fn add_assign(&mut self, rhs: State) {
        self.pos += rhs.pos;
        self.vel += rhs.vel;
    }
}

impl Mul<FloatType> for State {
    type Output = State;

    fn mul(self, rhs: FloatType) -> Self::Output {
        State {
            pos: self.pos * rhs,
            vel: self.vel * rhs,
        }
    }
}

impl MulAssign<FloatType> for State {
    fn mul_assign(&mut self, rhs: FloatType) {
        self.pos *= rhs;
        self.vel *= rhs;
    }
}
