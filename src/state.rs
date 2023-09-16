use crate::prelude::*;
use core::ops::{Add, AddAssign, Mul, MulAssign};

#[derive(Clone, Copy)]
pub struct State {
    pub pos: Vector3,
    pub vel: Vector3,
}

impl State {
    pub fn new(x0: Vector3, v0: Vector3) -> Self {
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
