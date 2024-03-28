//! # Single Operations Module
//!
//! The `single_ops` module implements unary operations for the `Perplex` struct.
//!
//! ## Traits from `num_traits` crate:
//! - `Neg`: Trait for negation, which returns the additive inverse of a perplex number.
//! - `Inv`: Trait for inversion, which provides the multiplicative inverse of a perplex number if it exists.
//! - `Pow`: Trait for both signed and unsigned integers, leveraging the exponentiation by squaring algorithm.
//!
//! ## Exponentiation Functions
//! The module defines methods for exponentiation:
//! - `powu`: Method for exponentiation with an unsigned integer exponent.
//! - `powi`: Method for exponentiation with a signed integer exponent, returning an `Option` to handle cases where the perplex number cannot be inverted.

use super::Perplex;
use num_traits::{Inv, Num, One, Pow};
use std::ops::Neg;

impl<T: Copy + Num + Neg<Output = T>> Neg for Perplex<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.t, -self.x)
    }
}
impl<T: Copy + Num + Neg<Output = T>> Inv for Perplex<T> {
    type Output = Option<Self>;
    #[inline]
    fn inv(self) -> Self::Output {
        self.try_inverse()
    }
}

impl<T: Copy + Num + Neg<Output = T>> Perplex<T> {
    /// Raises `self` to an unsigned integer power.
    #[inline]
    pub fn powu(&self, exp: u32) -> Self {
        Pow::pow(*self, exp)
    }

    /// Raises `self` to a signed integer power.
    #[inline]
    pub fn powi(&self, exp: i32) -> Option<Self> {
        Pow::pow(*self, exp)
    }
}
impl<T: Copy + Num> Pow<u32> for Perplex<T> {
    type Output = Perplex<T>;

    /// Performs exponentiation by squaring, an efficient algorithm for raising numbers to a power.
    /// This method is an iterative implementation of the algorithm described at [Exponentiation by Squaring](https://wikipedia.org/wiki/Exponentiation_by_squaring).
    ///
    /// # Arguments
    /// * `exp` - The exponent to raise the perplex number to.
    ///
    /// # Returns
    /// The result of raising the perplex number to the power of `exp`.
    #[inline]
    fn pow(self, mut exp: u32) -> Self::Output {
        // Initialize the result as the multiplicative identity, which is the result if the exponent is zero.
        let mut result = Perplex::one();
        if exp == 0 {
            return result;
        }
        // Set the base for exponentiation and iterate until the exponent is reduced to 1.
        let mut base = self;
        while exp > 1 {
            if exp % 2 == 1 {
                result = result * base;
            }
            exp /= 2;
            base = base * base;
        }
        result * base
    }
}

impl<T: Copy + Num + Neg<Output = T>> Pow<i32> for Perplex<T> {
    type Output = Option<Perplex<T>>;

    /// Performs exponentiation for both positive and negative integer exponents.
    /// For negative exponents, it calculates the multiplicative inverse before exponentiation.
    ///
    /// # Arguments
    /// * `exp` - The exponent to raise the perplex number to.
    ///
    /// # Returns
    /// An `Option` containing the result of raising the perplex number to the power of `exp`.
    /// Returns `None` if the perplex number cannot be inverted (i.e., it is light-like).
    #[inline]
    fn pow(self, exp: i32) -> Self::Output {
        // If the exponent is negative, calculate the multiplicative inverse first.
        if exp < 0 {
            // Use the wrapping_neg method to safely handle potential overflow.
            self.inv().map(|z| z.pow(exp.wrapping_neg() as u32))
        } else {
            // For non-negative exponents, delegate to the u32 implementation.
            Some(Pow::pow(self, exp as u32))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use num_traits::*;
    #[test]
    fn test_inv() {
        let z = Perplex::new(2.0, -1.0);
        let inv_result = z.inv();
        assert!(inv_result.is_some(), "2 - h is invertible!");
        assert_eq!(inv_result.unwrap(), Perplex::new(2.0 / 3.0, 1.0 / 3.0));
    }

    #[test]
    fn test_power_u32() {
        let z = Perplex::new(1.0, -1.0);
        assert_eq!(
            z.powu(0),
            Perplex::new(1.0, 0.0),
            "Power 0 yields neutral element of multiplication!"
        );
        assert_eq!(
            z * z,
            Perplex::new(2.0, -2.0),
            "Multiplication with itself!"
        );
        assert_eq!(
            z.powu(2),
            Perplex::new(2.0, -2.0),
            "Power 2 yields multiplication with itself!"
        );
        assert_eq!(
            z.powu(3),
            Perplex::new(4.0, -4.0),
            "Power 3 multiplication result!"
        );
        let z = Perplex::new(f64::PI(), -0.123);
        assert_eq!(z.powu(3), z * z * z, "Power 3 multiplication result!");
        assert_abs_diff_eq!(
            z.powu(8),
            z * z * z * z * z * z * z * z,
            epsilon = 0.0000001
        );
        assert_abs_diff_eq!(z.powu(7), z * z * z * z * z * z * z, epsilon = 0.0000001);
    }
    #[test]
    fn test_power_i32() {
        let z = Perplex::new(1.0, -1.0);
        assert!(z.powi(-2).is_none(), " 1 - h is not invertibe!");
        let z = Perplex::new(2.0, 1.0);
        let z_inv = z.try_inverse().unwrap();
        assert_eq!(
            z_inv.powi(2).unwrap(), // (2/3 - h / 3)^2
            Perplex::new(5.0 / 9.0, -4.0 / 9.0),
            "Multiplication of inverse with itself!"
        );
        assert_eq!(
            z.powi(-2).unwrap(),
            Perplex::new(5.0 / 9.0, -4.0 / 9.0),
            "Power -2 yields multiplication of inverse with itself!"
        );
        assert_eq!(
            z.powi(-3).unwrap(),
            Perplex::new(14.0 / 27.0, -13.0 / 27.0),
            "Power -3 multiplication result!"
        );
        let z = Perplex::new(f64::PI(), -0.123);
        let z_inv = z.try_inverse().unwrap();
        assert_eq!(
            z.powi(-3).unwrap(),
            z_inv * z_inv * z_inv,
            "Power -3 multiplication result!"
        );
        assert_abs_diff_eq!(
            z.powi(-8).unwrap(),
            z_inv * z_inv * z_inv * z_inv * z_inv * z_inv * z_inv * z_inv,
        );
        assert_abs_diff_eq!(
            z.powi(-7).unwrap(),
            z_inv * z_inv * z_inv * z_inv * z_inv * z_inv * z_inv,
        );
    }
}
