use crate::{prelude::*, state::State};

pub type OdeSolver = RK4;

pub struct RK4 {
    y: State,
    dt: FloatType,
    t: FloatType,
}

impl RK4 {
    pub fn new(y0: State, dt: FloatType) -> Self {
        RK4 { y: y0, dt, t: 0.0 }
    }

    pub fn step<F>(&mut self, func: F) -> (&State, FloatType)
    where
        F: Fn(State) -> State,
    {
        self.t += self.dt;

        let k1 = func(self.y);
        let k2 = func(self.y + k1 * (0.5 * self.dt));
        let k3 = func(self.y + k2 * (0.5 * self.dt));
        let k4 = func(self.y + k3 * self.dt);
        self.y += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * ((1.0 / 6.0) * self.dt);

        (&self.y, self.t)
    }
}
