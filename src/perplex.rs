//! # Perplex Module
//! This module defines the `Perplex` struct and provides common mathematical methods for it.
//!
//! ## Features
//! - Calculation of common distance metrics as well as the squared distance in the hyperbolic plane.
//! - Determination of the number's nature (time-like, space-like, or light-like) based on its squared distance. See Properties of the Perplex Numbers in [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643).
//! - `AbsDiffEq` trait from the `approx` crate.
//! - Constants and `FloatCore` traits from the `num_traits` crate.
//! - Hyperbolic exponential function as well as the natural logarithm as the inversion.
//! - Common trigonometric functions in the hyperbolic plane.

use approx::AbsDiffEq;
use num_traits::float::FloatCore;
use num_traits::{Float, Num, One, Zero};
use std::fmt;
use std::ops::Neg;

/// The `Perplex` struct is a representation of hyperbolic numbers, also known as split-complex numbers, which consist of two components: a real part (t) and a hyperbolic part (x). These components correspond to the time and space coordinates in Minkowski space-time, respectively. See Sec. 4.1 `Geometrical Representation of Hyperbolic Numbers` in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
/// The implementation is generic over a type `T`, which allows it to be used with different numeric types (i.e., `f32` or `f64`).
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Perplex<T> {
    /// The real part of the perplex number, representing time.
    pub t: T,
    /// The hyperbolic part of the perplex number, representing space.
    pub x: T,
}

impl<T> Perplex<T> {
    /// Create a new Perplex number
    #[inline]
    pub fn new(t: T, x: T) -> Self {
        Self { t, x }
    }
}

impl<T: Copy + Neg<Output = T> + PartialOrd + Num + fmt::Display> fmt::Display for Perplex<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (x, sign) = if self.x < T::zero() {
            (-self.x, "-")
        } else {
            (self.x, "+")
        };
        match f.precision() {
            Some(p) => write!(f, "{:.*} {sign} {:.*} h", p, self.t, p, x,),
            None => {
                let t_pretty = format!("{:.1$}", self.t, 2);
                let x_pretty = format!("{:.1$}", x, 2);
                write!(f, "{} {sign} {} h", t_pretty, x_pretty)
            }
        }
    }
}

impl<T: AbsDiffEq> AbsDiffEq for Perplex<T>
where
    T::Epsilon: Copy,
{
    type Epsilon = T::Epsilon;
    #[inline]
    fn default_epsilon() -> Self::Epsilon {
        T::default_epsilon()
    }
    #[inline]
    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        T::abs_diff_eq(&self.t, &other.t, epsilon) && T::abs_diff_eq(&self.x, &other.x, epsilon)
    }
}

impl<T: Copy + Num> Default for Perplex<T> {
    /// Defaults to the neutral element of multiplication.
    #[inline]
    fn default() -> Self {
        Self::new(T::one(), T::zero())
    }
}

impl<T: Copy + Num> From<T> for Perplex<T> {
    /// Conversion of a number `t` into a Perplex yields time-component `t` with zero space component.
    #[inline]
    fn from(t: T) -> Self {
        Self::new(t, T::zero())
    }
}

impl<T: Copy + Num> Perplex<T> {
    /// Returns hyperbolic unit.
    #[inline]
    pub fn h() -> Self {
        Self::new(T::zero(), T::one())
    }
    /// Returns the time component.
    #[inline]
    pub fn real(&self) -> T {
        self.t
    }
    /// Returns the space component.
    #[inline]
    pub fn hyperbolic(&self) -> T {
        self.x
    }
    /// Returns the squared distance D(z) in the hyperbolic plane.
    #[inline]
    pub fn squared_distance(&self) -> T {
        self.t * self.t - self.x * self.x
    }
    /// Multiplies `self` by the scalar `factor`.
    #[inline]
    pub fn scale(&self, factor: T) -> Self {
        Self::new(factor * self.t, factor * self.x)
    }
}
impl<T: Copy + Num + PartialOrd> Perplex<T> {
    /// Checks if the perplex number is time-like, i.e., the squared distance is positive.
    #[inline]
    pub fn is_time_like(&self) -> bool {
        self.squared_distance() > T::zero()
    }
    /// Checks if the perplex number is space-like, i.e., the squared distance is negative.
    #[inline]
    pub fn is_space_like(&self) -> bool {
        self.squared_distance() < T::zero()
    }
    /// Checks if the perplex number is light-like, i.e., the squared distance is zero.
    #[inline]
    pub fn is_light_like(&self) -> bool {
        self.squared_distance() == T::zero()
    }
}
impl<T: Copy + Num + Neg<Output = T>> Perplex<T> {
    /// Returns the hyperbolic conjugate.
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.t, -self.x)
    }
    /// Returns the multiplicative inverse `1/self`, if it exists, or `None` if not.
    #[inline]
    pub fn try_inverse(&self) -> Option<Self> {
        let squared_distance = self.squared_distance();
        if squared_distance == T::zero() {
            None
        } else {
            Some(Self::new(
                self.t / squared_distance,
                -self.x / squared_distance,
            ))
        }
    }
}

impl<T: Copy + Float> Perplex<T> {
    /// Returns the L1 norm `|t| + |x|` (Manhattan distance) from the origin in the cartesian coordinate plane, see Eq. 2.49 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
    #[inline]
    pub fn l1_norm(&self) -> T {
        self.t.abs() + self.x.abs()
    }
    /// Returns the L2 norm `|t^2| + |x^2|` (Euclidean distance) from the origin in the cartesian coordinate plane, see Eq. 2.50 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
    #[inline]
    pub fn l2_norm(&self) -> T {
        (self.t * self.t + self.x * self.x).sqrt()
    }
    /// Returns the maximum norm `||z||_∞ = max(|t|, |x|)` from the origin in the cartesian coordinate plane, see Eq. 2.51 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
    #[inline]
    pub fn max_norm(&self) -> T {
        self.t.abs().max(self.x.abs())
    }

    /// Returns the modulus of `self`.
    #[inline]
    pub fn modulus(self) -> T {
        let d_z = self.squared_distance();
        d_z.abs().sqrt()
    }
    /// Returns the norm (modulus) of `self`.
    #[inline]
    pub fn norm(self) -> T {
        self.modulus()
    }
    /// Returns the magnitude (modulus) of `self`.
    #[inline]
    pub fn magnitude(self) -> T {
        self.modulus()
    }

    /// Computes the hyperbolic exponential function for all sectors. Formula is extended to all sectors, see Sec 4.1.1 Hyperbolic Exponential Function and 7.4 The Elementary Functions of a Canonical Hyperbolic Variable in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn exp(self) -> Self {
        let k = self.klein().unwrap_or(Perplex::one());
        let Self { t, x } = k * self;
        let t_exp = t.exp();
        k * Self::new(t_exp * x.cosh(), t_exp * x.sinh())
    }
    /// Computes the inverse of the hyperbolic exponential function, i.e., the natural logarithm. Formula is extended to all sectors, see Sec. 7.4 The Elementary Functions of a Canonical Hyperbolic Variable in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn ln(self) -> Option<Self> {
        self.klein().map(|k| {
            let Self { t, x } = k * self;
            let squared_distance = t * t - x * x;
            let two = T::one() + T::one();
            let t_new = squared_distance.ln() / two;
            let x_new = (x / t).atanh();
            k * Self::new(t_new, x_new)
        })
    }

    /// Returns the logarithm of `self` with respect to an arbitrary base, if the natural logarithm of `self` exists, according to the formula `ln(self) / ln(base)`.
    #[inline]
    pub fn log(self, base: T) -> Option<Self> {
        self.ln().map(|z| z / base.ln())
    }

    /// Computes the square root of `self` if `self` lies in the right sector, or returns `None` if not. Formula is taken from Eq. 2.23 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
    #[inline]
    pub fn sqrt(self) -> Option<Self> {
        let t_x_add = self.t + self.x;
        let t_x_sub = self.t - self.x;
        if t_x_add >= T::zero() && t_x_sub >= T::zero() {
            let sqrt_add = t_x_add.sqrt();
            let sqrt_sub = t_x_sub.sqrt();
            let two = T::one() + T::one();
            let t = (sqrt_add + sqrt_sub) / two;
            let x = (sqrt_add - sqrt_sub) / two;
            Some(Perplex::new(t, x))
        } else {
            None
        }
    }

    /// Computes the sinus (circular trigonometric) of `self`. Formula is taken from Eq. 7.4.6 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn sin(self) -> Self {
        Self::new(self.t.sin() * self.x.cos(), self.t.cos() * self.x.sin())
    }
    /// Computes the cosinus (circular trigonometric) of `self`. Formula is taken from Eq. 7.4.6 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn cos(self) -> Self {
        Self::new(self.t.cos() * self.x.cos(), self.t.sin() * self.x.sin())
    }
    /// Computes the tangens (circular trigonometric) of `self` by the formula `sin(self) / cos(self)`. Returns `None` if `cos(self)` is light-like.
    #[inline]
    pub fn tan(self) -> Option<Self> {
        self.sin() / self.cos()
    }
    /// Computes the sinh (hyperbolic trigonometric) of `self`. Formula is taken from Eq. 7.4.5 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn sinh(self) -> Self {
        Self::new(self.t.sinh() * self.x.cosh(), self.t.cosh() * self.x.sinh())
    }
    /// Computes the cosh (hyperbolic trigonometric) of `self`. Formula is taken from Eq. 7.4.5 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn cosh(self) -> Self {
        Self::new(self.t.cosh() * self.x.cosh(), self.t.sinh() * self.x.sinh())
    }
    /// Computes the tanh (hyperbolic trigonometric) of `self` by the formula `sinh(self) / cosh(self)`. Returns `None` if `cosh(self)` is light-like.
    #[inline]
    pub fn tanh(self) -> Option<Self> {
        self.sinh() / self.cosh()
    }
}

impl<T: FloatCore> Perplex<T> {
    /// Checks if the given perplex number is NaN
    #[inline]
    pub fn is_nan(self) -> bool {
        self.t.is_nan() || self.x.is_nan()
    }

    /// Checks if the given perplex number is infinite
    #[inline]
    pub fn is_infinite(self) -> bool {
        !self.is_nan() && (self.t.is_infinite() || self.x.is_infinite())
    }

    /// Checks if the given perplex number is finite
    #[inline]
    pub fn is_finite(self) -> bool {
        self.t.is_finite() && self.x.is_finite()
    }

    /// Checks if the given perplex number is normal
    #[inline]
    pub fn is_normal(self) -> bool {
        self.t.is_normal() && self.x.is_normal()
    }
}

// constants
impl<T: Copy + Num> Zero for Perplex<T> {
    #[inline]
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.t.is_zero() && self.x.is_zero()
    }

    #[inline]
    fn set_zero(&mut self) {
        self.t.set_zero();
        self.x.set_zero();
    }
}

impl<T: Copy + Num> One for Perplex<T> {
    #[inline]
    fn one() -> Self {
        Self::new(One::one(), Zero::zero())
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.t.is_one() && self.x.is_zero()
    }

    #[inline]
    fn set_one(&mut self) {
        self.t.set_one();
        self.x.set_zero();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use num_traits::*;

    #[test]
    fn test_display() {
        let z = Perplex::new(1.1235, 1.10);
        assert_eq!(
            format!("{:.3}", z),
            String::from("1.123 + 1.100 h"),
            "Precision specifier produces 3 decimal places!"
        );
        assert_eq!(
            format!("{:.1}", z),
            String::from("1.1 + 1.1 h"),
            "Precision specifier produces 1 decimal place!"
        );
        let z = Perplex::new(2.0, -1.0);
        assert_eq!(z.to_string(), String::from("2.00 - 1.00 h"), "Negation sign is used for negative space component! Per default, fmt produces two decimal places!");
    }
    #[test]
    fn test_components() {
        let z = Perplex::new(1.1, 2.2);
        assert_eq!(z.real(), 1.1);
        assert_eq!(z.hyperbolic(), 2.2);
        assert_eq!(z.scale(2.0), Perplex::new(2.2, 4.4));
        assert_eq!(Perplex::from(2.0), Perplex::new(2.0, 0.0), "Converting a number t into a Perplex yields time-component t and zero space component!")
    }
    #[test]
    fn test_norm() {
        let z = Perplex::new(2.0, -1.0);
        assert!(z.is_time_like());
        assert_eq!(z.modulus(), f64::sqrt(3.0), "2 - h has a norm of √3");
        let z = Perplex::new(1.0, -1.0);
        assert!(z.is_light_like());
        assert_eq!(z.magnitude(), f64::zero(), "1 - h has a norm of zero");
        let z = Perplex::new(-1.0, 2.0);
        assert!(z.is_space_like());
        assert_eq!(z.norm(), f64::sqrt(3.0), "-1 + 2h has a norm of √3");
        assert_eq!(z.l1_norm(), 3.0, "-1 + 2h has a l1 norm of 3");
        assert_eq!(z.l2_norm(), f64::sqrt(5.0), "-1 + 2h has a l2 norm of √5");
        assert_eq!(z.max_norm(), 2.0, "-1 + 2h has a max norm of 2");
    }

    #[test]
    fn test_log() {
        let z = Perplex::new(2.0, 1.0);
        let z_ln = z.ln().unwrap();
        let z_log = z.log(2.0).unwrap();
        assert_eq!(z_log, z_ln / f64::ln(2.0));
    }
    #[test]
    fn test_logarithm_exponential() {
        let z = Perplex::new(2.0, 1.0); // Right-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for time-like 2 + h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for time-like -2 + h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(1.0, 2.0); // Up-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for space-like 1 + 2h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(1.0, -2.0); // Down-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for space-like 1 - 2h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);
    }
    #[test]
    fn test_exponential_logarithm() {
        let z = Perplex::new(2.0, 1.0); // Right-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(1.0, 2.0); // Up-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(1.0, -2.0); // Down-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
    }

    #[test]
    fn test_trigonometric() {
        let pi = f64::PI();
        let z = Perplex::new(pi, pi / 2.0).sin();
        assert_abs_diff_eq!(z, Perplex::new(0.0, -1.0));
        assert!(
            z.tan().is_some(),
            "Tangens of z should be defined since cos(z) is not light-like!"
        );
        let zero: Perplex<f64> = Perplex::zero();
        assert_abs_diff_eq!(zero.sinh(), zero);
        let z = Perplex::new(1.0, 0.0);
        let expected_tanh = Perplex::new(z.t.tanh(), 0.0);
        assert_eq!(
            z.tanh(),
            Some(expected_tanh),
            "Tanh of z should be defined since cosh(z) is not light-like!"
        );
    }

    #[test]
    fn test_sqrt() {
        // Test sqrt for a Perplex number in the Right sector (t > |x|)
        let z_right = Perplex::new(2.0, 1.0);
        assert!(
            z_right.sqrt().is_some(),
            "Sqrt should be defined for Perplex numbers in the Right sector."
        );
        // The expected result should be a Perplex number whose square equals z_right
        if let Some(sqrt_z) = z_right.sqrt() {
            assert_abs_diff_eq!(sqrt_z.powu(2), z_right, epsilon = 1e-10);
        }
        // Test sqrt for a Perplex number in the Left sector (t < -|x|)
        let z_left = Perplex::new(-2.0, 1.0);
        assert!(
            z_left.sqrt().is_none(),
            "Sqrt should not be defined for Perplex numbers in the Left sector."
        );
    }

    #[test]
    fn test_core() {
        let z = Perplex::new(1.0, 2.0);
        assert!(
            z.is_finite(),
            "Perplex number with finite components is finite!"
        );
        assert!(
            z.is_normal(),
            "Perplex number with finite components is normal!"
        );
        assert!(
            !z.is_infinite(),
            "Perplex number with finite components is not infinite!"
        );
        assert!(
            !z.is_nan(),
            "Perplex number with finite components is not NAN!"
        );
        let z = Perplex::new(f64::NAN, 1.0);
        assert!(z.is_nan(), "Perplex number with a NaN component is NAN!")
    }
    #[test]
    fn test_const() {
        let mut z = Perplex::new(0.0, 2.0);
        assert!(!z.is_one(), "Perplex number a zero component is not one!");
        assert!(
            !z.is_zero(),
            "Perplex number with non-zero component is not zero!"
        );
        z.set_one();
        assert_eq!(Perplex::one(), z, "Perplex number set to one equals one!");
        assert!(z.is_one(), "Perplex number set to one is one!");
        z.set_zero();
        assert_eq!(
            Perplex::zero(),
            z,
            "Perplex number set to zero equals zero!"
        );
        assert!(z.is_zero(), "Perplex number set to zero is zero!");

        assert!(
            !z.is_infinite(),
            "Perplex number with finite components is not infinite!"
        );
        assert!(
            !z.is_nan(),
            "Perplex number with finite components is not NAN!"
        );
        let z = Perplex::new(f64::NAN, 1.0);
        assert!(z.is_nan(), "Perplex number with a NaN component is NAN!")
    }
}
