use crate::state::State;

pub struct RK4 {
    y: State,
    dt: f64,
}

impl RK4 {
    pub fn new(y0: State, dt: f64) -> Self {
        RK4 { y: y0, dt }
    }

    pub fn step<F>(&mut self, func: F)
    where
        F: Fn(State) -> State,
    {
        let k1 = func(self.y);
        let k2 = func(self.y + k1 * (0.5 * self.dt));
        let k3 = func(self.y + k2 * (0.5 * self.dt));
        let k4 = func(self.y + k3 * self.dt);
        self.y += (k1 +  k2 * 2.0 + k3 * 2.0 + k4) * ((1.0 / 6.0) * self.dt);
    }

    pub fn get_state(&self) -> &State {
        &self.y
    }
}