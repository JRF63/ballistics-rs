use crate::{data::cubic_hermite_interpolation, prelude::*, state::State};
use core::{
    cmp::Ordering,
    ops::{Add, Mul},
};

pub type OdeSolver = RK4;

pub struct RK4 {
    y: State,
    yp: State,
    t: FloatType,
    y_prev: State,
    yp_prev: State,
    t_prev: FloatType,
    dt: FloatType,
}

impl RK4 {
    pub fn new<F>(y0: State, func: F, dt: FloatType) -> Self
    where
        F: Fn(State) -> State,
    {
        let yp = func(y0);
        Self {
            y: y0,
            yp,
            t: 0.0,
            y_prev: y0,
            yp_prev: yp,
            t_prev: 0.0,
            dt,
        }
    }

    pub fn step<F>(&mut self, func: F)
    where
        F: Fn(State) -> State,
    {
        self.y_prev = self.y;
        self.yp_prev = self.yp;
        self.t_prev = self.t;

        let k1 = self.yp;
        let k2 = func(self.y + k1 * (0.5 * self.dt));
        let k3 = func(self.y + k2 * (0.5 * self.dt));
        let k4 = func(self.y + k3 * self.dt);

        self.y += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * ((1.0 / 6.0) * self.dt);
        self.yp = func(self.y);
        self.t += self.dt;
    }

    pub fn current_state(&self) -> (&State, FloatType) {
        (&self.y, self.t)
    }

    pub fn find_state_at_event<E, F, G>(&self, extractor: F, event: G) -> Option<(State, FloatType)>
    where
        F: Fn(&State) -> E + Copy,
        G: Fn(E) -> FloatType + Copy,
        E: Add<Output = E> + Mul<FloatType, Output = E> + Copy,
    {
        const MAX_ITERS: u32 = 64;

        let y0 = extractor(&self.y_prev);
        let y1 = extractor(&self.y);
        let f0 = extractor(&self.yp_prev);
        let f1 = extractor(&self.yp);

        let event_prev = event(y0);
        let event_now = event(y1);

        // If the target event is reached at the boundaries, return them immediately
        if event_prev == 0.0 {
            return Some((self.y_prev, self.t_prev));
        }
        if event_now == 0.0 {
            return Some((self.y, self.t));
        }

        // One end should be positive while the other is negative for the zero crossing
        if event_prev * event_now < 0.0 {
            let mut low = 0.0;
            let mut high = 1.0;
            let mut theta = (low + high) / 2.0;

            let low_is_negative = event_prev.is_sign_negative();

            for _ in 0..MAX_ITERS {
                let root = cubic_hermite_interpolation(y0, y1, f0, f1, self.dt, theta);
                match event(root).total_cmp(&0.0) {
                    Ordering::Less => {
                        if low_is_negative {
                            low = theta;
                        } else {
                            high = theta;
                        }
                    }
                    Ordering::Greater => {
                        if low_is_negative {
                            high = theta;
                        } else {
                            low = theta;
                        }
                    }
                    Ordering::Equal => break,
                }

                let theta_old = theta;
                theta = (low + high) / 2.0;
                if theta == theta_old {
                    break;
                }
            }

            let y_interpolated = cubic_hermite_interpolation(
                self.y_prev,
                self.y,
                self.yp_prev,
                self.yp,
                self.dt,
                theta,
            );
            let t_interpolated = self.t_prev + self.dt * theta;
            Some((y_interpolated, t_interpolated))
        } else {
            None
        }
    }
}
