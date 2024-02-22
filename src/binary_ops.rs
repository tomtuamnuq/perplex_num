use super::Perplex;
use num_traits::{Num, NumAssign};
use std::ops::{Add, Div, Mul, Sub};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

// Properties of the Perplex Numbers in Fundamental Theorems of Algebra for the Perplexes:
// binary between Perplex and Perplex
impl<T: Copy + Num> Add for Perplex<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.t + rhs.t, self.x + rhs.x)
    }
}
impl<T: Copy + NumAssign> AddAssign for Perplex<T> {
    fn add_assign(&mut self, rhs: Self) {
        self.t += rhs.t;
        self.x += rhs.x;
    }
}

impl<T: Copy + Num> Sub for Perplex<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.t - rhs.t, self.x - rhs.x)
    }
}
impl<T: Copy + NumAssign> SubAssign for Perplex<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self.t -= rhs.t;
        self.x -= rhs.x;
    }
}

impl<T: Copy + Num> Mul for Perplex<T> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.t * rhs.t + self.x * rhs.x,
            rhs.t * self.x + self.t * rhs.x,
        )
    }
}
impl<T: Copy + NumAssign> MulAssign for Perplex<T> {
    fn mul_assign(&mut self, rhs: Self) {
        let t = self.t;
        self.t *= rhs.t;
        self.t += self.x * rhs.x;
        self.x *= rhs.t;
        self.x += t * rhs.x;
    }
}

impl<T: Copy + Num> Div for Perplex<T> {
    type Output = Option<Self>;
    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let Self { t: t2, x: x2 } = rhs;
        let norm_squared_2 = t2 * t2 - x2 * x2;
        if norm_squared_2 == T::zero() {
            // light-like
            None
        } else {
            let Self { t: t1, x: x1 } = self;
            let t_new = (t1 * t2 - x1 * x2) / norm_squared_2;
            let x_new = (t2 * x1 - t1 * x2) / norm_squared_2;
            Some(Self::new(t_new, x_new))
        }
    }
}
impl<T: Copy + NumAssign> DivAssign for Perplex<T> {
    fn div_assign(&mut self, rhs: Self) {
        let Self { t: t2, x: x2 } = rhs;
        let norm_squared_2 = t2 * t2 - x2 * x2;
        if norm_squared_2 == T::zero() {
            // light-like
            return;
        }
        let t = self.t;
        self.t *= t2;
        self.t -= self.x * x2;
        self.t /= norm_squared_2;
        self.x *= t2;
        self.x -= t * x2;
        self.x /= norm_squared_2;
    }
}

// TODO binary ops between references ? use copy semantics

// binary between Perplex and T
impl<T: Copy + Num> Add<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn add(self, other: T) -> Self::Output {
        Self::Output::new(self.t + other, self.x)
    }
}
impl<T: Copy + NumAssign> AddAssign<T> for Perplex<T> {
    fn add_assign(&mut self, rhs: T) {
        self.t += rhs;
    }
}

impl<T: Copy + Num> Sub<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        Self::Output::new(self.t - rhs, self.x)
    }
}
impl<T: Copy + NumAssign> SubAssign<T> for Perplex<T> {
    fn sub_assign(&mut self, rhs: T) {
        self.t -= rhs;
    }
}

impl<T: Copy + Num> Mul<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Self::Output::new(self.t * rhs, self.x * rhs)
    }
}
impl<T: Copy + NumAssign> MulAssign<T> for Perplex<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.t *= rhs;
        self.x *= rhs;
    }
}

impl<T: Copy + Num> Div<T> for Perplex<T> {
    type Output = Self;
    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Self::Output::new(self.t / rhs, self.x / rhs)
    }
}
impl<T: Copy + NumAssign> DivAssign<T> for Perplex<T> {
    fn div_assign(&mut self, rhs: T) {
        self.t /= rhs;
        self.x /= rhs;
    }
}
