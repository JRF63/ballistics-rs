use std::ops::{Add, AddAssign, Mul, MulAssign};

use glam::f64::DVec3;

#[derive(Clone, Copy)]
pub struct State(pub [DVec3; 2]);

impl Add<State> for State {
    type Output = State;

    fn add(self, rhs: State) -> Self::Output {
        State([self.0[0] + rhs.0[0], self.0[1] + rhs.0[1]])
    }
}

impl AddAssign<State> for State {
    fn add_assign(&mut self, rhs: State) {
        self.0[0] += rhs.0[0];
        self.0[1] += rhs.0[1];
    }
}

impl Mul<f64> for State {
    type Output = State;

    fn mul(self, rhs: f64) -> Self::Output {
        State([self.0[0] * rhs, self.0[1] * rhs])
    }
}

impl MulAssign<f64> for State {
    fn mul_assign(&mut self, rhs: f64) {
        self.0[0] *= rhs;
        self.0[1] *= rhs;
    }
}