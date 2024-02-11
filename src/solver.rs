use crate::{data::cubic_hermite_interpolation, prelude::*, state::State};
use core::ops::{Add, Mul};

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

    pub fn current_state(&self) -> (&State, &State, FloatType) {
        (&self.y, &self.yp, self.t)
    }

    pub fn find_state_at_event<E, F, G>(&self, extractor: F, event: G) -> Option<(State, FloatType)>
    where
        F: Fn(&State) -> E + Copy,
        G: Fn(E) -> FloatType + Copy,
        E: Add<Output = E> + Mul<FloatType, Output = E> + Copy,
    {
        const MAX_ITERS: i32 = 64;

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
            let func = |theta: FloatType| -> FloatType {
                event(cubic_hermite_interpolation(y0, y1, f0, f1, self.dt, theta))
            };

            let theta = brentq(func, 0.0, 1.0, 2e-12, 4.0 * FloatType::EPSILON, MAX_ITERS);

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

// Port of
// https://github.com/scipy/scipy/blob/ae3ca70a72ef04547fbc28501009246cef5ee6c8/scipy/optimize/Zeros/brentq.c#L37
fn brentq<F>(
    func: F,
    xa: FloatType,
    xb: FloatType,
    xtol: FloatType,
    rtol: FloatType,
    iter: i32,
) -> FloatType
where
    F: Fn(FloatType) -> FloatType + Copy,
{
    let mut xpre = xa;
    let mut xcur = xb;

    let mut xblk = 0.0;
    let mut fblk = 0.0;
    let mut spre = 0.0;
    let mut scur = 0.0;

    let mut fpre = func(xpre);
    let mut fcur = func(xcur);

    if fpre == 0.0 {
        return xpre;
    }
    if fcur == 0.0 {
        return xcur;
    }
    if fpre.signum() == fcur.signum() {
        return 0.0;
    }

    for _ in 0..iter {
        if fpre != 0.0 && fcur != 0.0 && fpre.signum() != fcur.signum() {
            xblk = xpre;
            fblk = fpre;
            scur = xcur - xpre;
            spre = scur;
        }
        if fblk.abs() < fcur.abs() {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        /* the tolerance is 2*delta */
        let delta = (xtol + rtol * xcur.abs()) / 2.0;
        let sbis = (xblk - xcur) / 2.0;

        if fcur == 0.0 || sbis.abs() < delta {
            return xcur;
        }

        if spre.abs() > delta && fcur.abs() < fpre.abs() {
            let stry = if xpre == xblk {
                /* interpolate */
                -fcur * (xcur - xpre) / (fcur - fpre)
            } else {
                /* extrapolate */
                let dpre = (fpre - fcur) / (xpre - xcur);
                let dblk = (fblk - fcur) / (xblk - xcur);
                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            };

            if 2.0 * stry.abs() < FloatType::min(spre.abs(), 3.0 * sbis.abs() - delta) {
                /* good short step */
                spre = scur;
                scur = stry;
            } else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }
        } else {
            /* bisect */
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;
        if scur.abs() > delta {
            xcur += scur;
        } else {
            xcur += if sbis > 0.0 { delta } else { -delta };
        }

        fcur = func(xcur);
    }

    xcur
}
