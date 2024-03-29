//! # Binary Operations Module
//!
//! This module is dedicated to implementing traits from both the standard library and the `num_traits` crate for the `Perplex` struct. It provides the foundational binary operations necessary for working with two perplex numbers. See `Properties of the Perplex Numbers` in [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
//!
//! The module includes implementations for basic arithmetic operations between `Perplex` structs, such as:
//! - `Add`: Trait for the addition operator.
//! - `Sub`: Trait for the subtraction operator.
//! - `Mul`: Trait for the multiplication operator.
//! - `Div`: Trait for the division operator.
//! - Tertiary operation `MulAdd` from the `num_traits` crate.
//!
//! Additionally, it supports assignment variants of these operations for mutable references of `Perplex` structs, which are:
//! - `AddAssign`: Trait for addition assignment.
//! - `SubAssign`: Trait for subtraction assignment.
//! - `MulAssign`: Trait for multiplication assignment.
//! - `DivAssign`: Trait for division assignment.
//! - Tertiary operation `MulAddAssign` from the `num_traits` crate.
//!
//! The module also includes implementations for interactions between `Perplex` structs and the generic floating point type (`f32` or `f64`).

use super::Perplex;
use num_traits::{MulAdd, MulAddAssign, Num, NumAssign};
use std::ops::{Add, Div, Mul, Sub};
use std::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

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
    /// Divides `self` by `rhs`. Division by a light-like number yields `None`, otherwise `Some(self / rhs)`.
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
    /// Divides `self` by `rhs` in place. Division by a light-like number yields a Perplex number with NaN components.
    fn div_assign(&mut self, rhs: Self) {
        let Self { t: t2, x: x2 } = rhs;
        let norm_squared_2 = t2 * t2 - x2 * x2;
        let t = self.t;
        self.t *= t2;
        self.t -= self.x * x2;
        self.t /= norm_squared_2;
        self.x *= t2;
        self.x -= t * x2;
        self.x /= norm_squared_2;
    }
}

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

// tertiary ops between three Perplex
impl<T: Copy + Num + MulAdd<Output = T>> MulAdd<Perplex<T>> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn mul_add(self, other: Perplex<T>, add: Perplex<T>) -> Self {
        let t = self.t * other.t + self.x * other.x + add.t;
        let x = other.t * self.x + self.t * other.x + add.x;
        Self::new(t, x)
    }
}
impl<T: Copy + NumAssign + MulAddAssign> MulAddAssign for Perplex<T> {
    fn mul_add_assign(&mut self, other: Self, add: Self) {
        let t = self.t;
        self.t *= other.t;
        self.t += self.x * other.x + add.t;
        self.x *= other.t;
        self.x += t * other.x + add.x;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::*;
    #[test]
    fn test_add() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 + one + zero,
            Perplex::new(2.0, 2.0),
            "Componentwise addition!"
        );
        assert_eq!(
            z1 + z1.conj(),
            Perplex::new(2.0, 0.0),
            "Addition of conjugate zeros the hyperbolic part!"
        );
        let mut z2 = Perplex::new(-3.0, 2.0);
        let z12 = z1 + z2;
        z2 += z1;
        assert_eq!(z12, z2, "AddAssign yields same result as Add!");
    }
    #[test]
    fn test_sub() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 - one - zero,
            Perplex::new(0.0, 2.0),
            "Componentwise subtraction!"
        );
        assert_eq!(
            z1 - z1.conj(),
            Perplex::new(0.0, 4.0),
            "Subtraction of conjugate doubles the hyperbolic part!"
        );
        let mut z2 = Perplex::new(-3.0, 2.0);
        let z12 = z2 - z1;
        z2 -= z1;
        assert_eq!(z12, z2, "SubAssign yields same result as Sub!");
    }
    #[test]
    fn test_mul() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 * one,
            z1,
            "Neutral element of multiplication yields same element!"
        );
        assert_eq!(z1 * zero, zero, "Neutral element of addition yields zero!");
        let mut z2 = Perplex::new(-1.0, 2.0);
        let z12 = z1 * z2;
        z2 *= z1;
        assert_eq!(z2, Perplex::new(3.0, 0.0), "Multiplication formula!");
        assert_eq!(z12, z2, "MulAssign yields same result as Mul!");
    }
    #[test]
    fn test_div() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            (z1 / one).unwrap(),
            z1,
            "Division of neutral element of multiplication yields same element!"
        );
        assert!(
            (z1 / zero).is_none(),
            "Division of neutral element of addition yields none!"
        );
        let z2 = Perplex::new(-1.0, 2.0);
        let mut z12 = z1 * z2;
        let div_result = z12 / z2;
        assert!(
            div_result.is_some(),
            "Division of product by multiplier is valid!"
        );
        assert_eq!(
            div_result.unwrap(),
            z1,
            "Division of product by multiplier gives multiplicand."
        );
        z12 /= z2;
        assert_eq!(z12, z1, "DivAssign yields same result as Div!");

        let z2 = Perplex::new(-1.0, 1.0);
        let mut z12 = z1 * z2;
        assert_eq!(z12, Perplex::new(1.0, -1.0), "Multiplication formula!");
        assert!(z2.is_light_like(), "-1 + j is light-like!");
        assert!(
            (z12 / z2).is_none(),
            "Division is not defined for light-like numbers!"
        );
        z12 /= z2;
        assert!(
            z12.t.is_nan() && z12.x.is_nan(),
            "DivAssign for light-like number yields NaN!"
        );
    }
    #[test]
    fn test_scalar() {
        let z1 = Perplex::new(1.0, 2.0);
        assert_eq!(
            z1 + 2.0,
            Perplex::new(3.0, 2.0),
            "Addition of scalar only on time component!"
        );
        assert_eq!(
            z1 - 2.0,
            Perplex::new(-1.0, 2.0),
            "Subtraction of scalar only on time component!"
        );
        assert_eq!(
            z1 * 2.0,
            Perplex::new(2.0, 4.0),
            "Componentwise scalar multiplication!"
        );
        assert_eq!(
            z1 / 2.0,
            Perplex::new(0.5, 1.0),
            "Componentwise scalar division!"
        );
    }
    #[test]
    fn test_scalar_assign() {
        let mut z1 = Perplex::new(1.0, 2.0);
        z1 += 2.0;
        assert_eq!(
            z1,
            Perplex::new(3.0, 2.0),
            "AddAssign of scalar only on time component!"
        );
        z1 -= 2.0;
        assert_eq!(
            z1,
            Perplex::new(1.0, 2.0),
            "SubAssign of scalar only on time component!"
        );
        z1 *= 2.0;
        assert_eq!(
            z1,
            Perplex::new(2.0, 4.0),
            "MulAssign componentwise scalar multiplication!"
        );
        z1 /= 2.0;
        assert_eq!(
            z1,
            Perplex::new(1.0, 2.0),
            "DivAssign componentwise scalar division!"
        );
    }
    #[test]
    fn test_mul_add() {
        let mut z1 = Perplex::new(1.0, 2.0);
        let z_mul = Perplex::new(-1.0, 2.0);
        let z_add = Perplex::new(-2.0, 1.0);
        let z = z1.mul_add(z_mul, z_add);
        z1.mul_add_assign(z_mul, z_add);
        assert_eq!(
            z,
            Perplex::new(1.0, 1.0),
            "Multiplication formula and addition!"
        );
        assert_eq!(z, z1, "MulAddAssign yields same result as MulAdd!");
    }
}
