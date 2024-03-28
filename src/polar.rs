//! # Hyperbolic Polar Module
//!
//! This module provides the functionality to work with perplex numbers in polar form, which is particularly useful in the context of hyperbolic geometry.
//! It includes methods for converting between the standard `Perplex` representation and the `HyperbolicPolar` form, as well as operations like exponentiation within the hyperbolic plane.
//! The hyperbolic polar form encodes a perplex number `z` as a triple of two real numbers `rho` and `theta`, as well as one out of four perplex numbers `klein`, such that `z= klein rho (cosh(theta) + h sinh(theta))`.
//! `Klein` is defined by the sector of the hyperbolic plane in which the perplex number is in. Formulas are taken from Tab. 1 and Appendix B in [Hyperbolic trigonometry in two-dimensional space-time geometry](https://doi.org/10.1393/ncb/i2003-10012-9).
//!
//! ## Usage
//!
//! Here is an example of how to use the `HyperbolicPolar` struct to convert a `Perplex` number
//! into its polar form and back:
//!
//! ```
//! use num_traits::Pow;
//! use perplex_num::{HyperbolicSector, HyperbolicPolar, Perplex};
//! let z = Perplex { t: 1.0, x: 0.5 };
//! // Convert the Perplex number to HyperbolicPolar form
//! let polar_form: HyperbolicPolar<f64> = z.into();
//! // Perform operations in polar form...
//! // For example, raise to a power
//! let polar_powered = polar_form.pow(2);
//! // Convert back to Perplex form
//! let z_powered: Perplex<f64> = polar_powered.into();
//! approx::assert_abs_diff_eq!(z_powered, Perplex { t: 1.25, x: 1.0 }, epsilon=0.0000000001);
//! ```

use super::Perplex;
use num_traits::{Float, Num, One, Pow};

/// Represents the sector of the hyperbolic plane a perplex number is in.
///
/// The hyperbolic plane is divided into four sectors by the intersection of two diagonals,
/// where `t = x` and `t = -x`. This enum also includes the `Diagonal` variant to represent
/// light-like perplex numbers where the time and space components are equal in magnitude.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub enum HyperbolicSector<T> {
    /// The sector where the time component is greater than the space component in absolute value.
    #[default]
    Right,
    /// The sector where the space component is greater than the time component in absolute value.
    Up,
    /// The mirror image of the Right sector, where the negative time component is greater in magnitude.
    Left,
    /// The mirror image of the Up sector, where the negative space component is greater in magnitude.
    Down,
    /// Represents a light-like perplex number on the diagonal where `t` and `x` are equal.
    /// The value `T` encodes which diagonal line is used based on its sign.
    Diagonal(T),
}

impl<T: Copy + Float> From<Perplex<T>> for HyperbolicSector<T> {
    /// Converts a perplex number into its corresponding hyperbolic sector.
    ///
    /// Light-like numbers are converted to the `Diagonal` variant, while others are
    /// categorized based on the magnitude and sign of their time and space components.
    #[inline]
    fn from(z: Perplex<T>) -> Self {
        let Perplex { t, x } = z;
        let (t_abs, x_abs) = (t.abs(), x.abs());
        if t_abs == x_abs {
            Self::Diagonal(t)
        } else if t_abs > x_abs {
            if t > T::zero() {
                Self::Right
            } else {
                Self::Left
            }
        } else if x > T::zero() {
            Self::Up
        } else {
            Self::Down
        }
    }
}

/// Represents a perplex number in hyperbolic polar form.
///
/// This struct is used to convert a perplex number to and from hyperbolic polar form,
/// which is useful for operations that are more naturally expressed in this form.
/// The conversion formulas are based on hyperbolic trigonometry principles.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct HyperbolicPolar<T> {
    /// The modulus of the perplex number, representing the hyperbolic radius.
    pub rho: T,
    /// The argument of the perplex number, representing the hyperbolic angle.
    pub theta: T,
    /// The sector of the hyperbolic plane the number is in.
    pub sector: HyperbolicSector<T>,
}

impl<T: Copy + Num> Default for HyperbolicPolar<T> {
    /// Provides a default value for `HyperbolicPolar`, which is in the Right sector
    /// with a modulus of one and an angle of zero, i.e., the neutral element of multiplication `Perplex::default()`.
    #[inline]
    fn default() -> Self {
        Self {
            rho: T::one(),
            theta: T::zero(),
            sector: HyperbolicSector::Right,
        }
    }
}

impl<T: Copy + Float> From<Perplex<T>> for HyperbolicPolar<T> {
    /// Converts a perplex number to hyperbolic polar form.
    ///
    /// The conversion takes into account the sector of the hyperbolic plane the number
    /// is in and uses the appropriate hyperbolic trigonometric functions.
    #[inline]
    fn from(z: Perplex<T>) -> Self {
        Self {
            rho: z.norm(),
            theta: z.arg(),
            sector: HyperbolicSector::from(z),
        }
    }
}

impl<T: Copy + Float> From<HyperbolicPolar<T>> for Perplex<T> {
    /// Converts a hyperbolic polar representation back to a perplex number.
    ///
    /// This method applies the inverse of the hyperbolic polar conversion, reconstructing
    /// the original perplex number from its modulus, argument, and sector.
    #[inline]
    fn from(polar: HyperbolicPolar<T>) -> Self {
        let HyperbolicPolar { rho, theta, sector } = polar;
        match sector {
            HyperbolicSector::Right => Self::new(rho * theta.cosh(), rho * theta.sinh()),
            HyperbolicSector::Up => Self::new(rho * theta.sinh(), rho * theta.cosh()),
            HyperbolicSector::Left => Self::new(-rho * theta.cosh(), -rho * theta.sinh()),
            HyperbolicSector::Down => Self::new(-rho * theta.sinh(), -rho * theta.cosh()),
            HyperbolicSector::Diagonal(t) => {
                if theta == T::infinity() {
                    Self::new(t, t)
                } else {
                    // theta == T::neg_infinity()
                    Self::new(t, -t)
                }
            }
        }
    }
}

impl<T: Copy + Float> Perplex<T> {
    /// Creates a new `Perplex` number `z`  with a given hyperbolic angle `theta` such that `z= exp(h theta)`.
    ///
    /// It is used to create a `Perplex` number with a given phase, using hyperbolic cosine and sine.
    #[inline]
    pub fn cis(theta: T) -> Self {
        Self::new(theta.cosh(), theta.sinh())
    }

    /// Calculates the hyperbolic argument of `self`.
    ///
    /// The argument is the angle in the hyperbolic plane from the positive time axis to the line
    /// connecting the origin to `self`. It is defined by a piecewise function, with special cases
    /// for light-like perplex numbers, whereby lines x=t and x=-t are mapped to ∞ and -∞, respectively, according to Sec. 4.1 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
    /// The formula is taken from Eq. 4.1.6 in Sec 4.1.1 `Hyperbolic Exponential Function and Hyperbolic Polar Transformation` in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn arg(self) -> T {
        let Self { t, x } = self;
        let (t_abs, x_abs) = (t.abs(), x.abs());
        if t_abs == x_abs {
            // self.is_light_like()
            if t == x {
                // line x = t
                T::infinity()
            } else {
                // line x = -t
                T::neg_infinity()
            }
        } else if t_abs > x_abs {
            (x / t).atanh()
        } else {
            (t / x).atanh()
        }
    }

    /// Calculate the Klein index of `self` for space- or time-like numbers. Returns `None` for light-like numbers.
    ///
    /// The Klein index is determined by the sector of the hyperbolic plane in which `self` resides.
    /// Formula is taken from Tab. 1 and Appendix B in [Hyperbolic trigonometry in two-dimensional space-time geometry](https://doi.org/10.1393/ncb/i2003-10012-9).
    #[inline]
    pub fn klein(self) -> Option<Self> {
        let Self { t, x } = self;
        let (t_abs, x_abs) = (t.abs(), x.abs());
        if t_abs == x_abs {
            // light-like
            None
        } else if t_abs > x_abs {
            if t > T::zero() {
                // Right-Sector
                Some(Self::one())
            } else {
                // Left-Sector
                Some(-Self::one())
            }
        } else if x > T::zero() {
            // Up-Sector
            Some(Self::h())
        } else {
            // Down-Sector
            Some(-Self::h())
        }
    }

    /// Retrieves the hyperbolic sector of the perplex number.
    ///
    /// # Examples
    ///
    /// ```
    /// use perplex_num::{Perplex, HyperbolicSector};
    ///
    /// let perplex = Perplex::new(1.0, 0.5);
    /// assert_eq!(perplex.sector(), HyperbolicSector::Right);
    /// ```
    #[inline]
    pub fn sector(&self) -> HyperbolicSector<T> {
        (*self).into()
    }

    /// Retrieves the hyperbolic polar form from a perplex number.
    ///
    /// # Examples
    ///
    /// ```
    /// use perplex_num::{Perplex, HyperbolicPolar};
    ///
    /// let perplex = Perplex::new(1.0, 0.5);
    /// let polar = perplex.polar();
    /// assert_eq!(polar.rho, perplex.norm());
    /// assert_eq!(polar.theta, perplex.arg());
    /// ```
    #[inline]
    pub fn polar(&self) -> HyperbolicPolar<T> {
        (*self).into()
    }
}

impl<T: Copy + Float> Pow<u32> for HyperbolicPolar<T> {
    /// Raises `self` to the power of unsigned `exp`.
    ///
    /// This method is based on an extended version of Formula 4.6 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586), ensuring consistency across the plane.
    type Output = Self;
    #[inline]
    fn pow(self, exp: u32) -> Self::Output {
        match exp {
            0 => Self::default(),
            1 => self,
            _ => {
                let n = exp as i32;
                let Self { rho, theta, sector } = self;
                if let HyperbolicSector::Diagonal(t) = sector {
                    let t_new = t * (t + t).powi(n - 1); // t^n * 2^{n-1}
                    HyperbolicPolar {
                        rho,
                        theta,
                        sector: HyperbolicSector::Diagonal(t_new),
                    }
                } else {
                    let new_sector = if n % 2 == 0 {
                        HyperbolicSector::Right // since -1^2 = 1 and h^2=1
                    } else {
                        sector
                    };
                    HyperbolicPolar {
                        rho: rho.powi(n), // Formula 4.6
                        theta: T::from(n).unwrap() * theta,
                        sector: new_sector,
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use num_traits::*;
    #[test]
    fn test_polar() {
        assert_abs_diff_eq!(
            Perplex::<f64>::default(),
            Perplex::from(HyperbolicPolar::default()),
            epsilon = 0.0001
        );

        let z = Perplex::new(1.0, 1.0); // Diagonal x=t
        assert!(z.is_light_like(), "1 + h is light-like!");
        assert_eq!(z.arg(), f64::infinity(), "Argument of 1 + h is infinity!");
        assert!(
            z.klein().is_none(),
            "Klein is not defined for light-like numbers!"
        );
        assert_eq!(
            HyperbolicPolar::from(z),
            HyperbolicPolar {
                rho: 0.0,
                theta: f64::infinity(),
                sector: HyperbolicSector::Diagonal(1.0)
            },
            "Polar form of 1 + h!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::new(1.0, -1.0); // Diagonal x=-t
        assert!(z.is_light_like(), "1 - h is light-like!");
        assert_eq!(
            z.arg(),
            f64::neg_infinity(),
            "Argument of 1 - h is  negative infinity!"
        );
        assert!(
            z.klein().is_none(),
            "Klein is not defined for light-like numbers!"
        );
        assert_eq!(
            HyperbolicPolar::from(z),
            HyperbolicPolar {
                rho: 0.0,
                theta: f64::neg_infinity(),
                sector: HyperbolicSector::Diagonal(1.0)
            },
            "Polar form of 1 - h!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::new(2.0, 1.0); // Right-Sector
        assert!(z.is_time_like(), "2 + h is time-like!");
        assert_ne!(z.arg(), 0.0, "2 + h has a non-zero argument!");
        assert_eq!(
            z.klein().unwrap(),
            Perplex::new(1.0, 0.0),
            "2 + h is in the right-sector!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        assert!(z.is_time_like(), "-2 + h is time-like!");
        assert_ne!(z.arg(), 0.0, "-2 + h has a non-zero argument!");
        assert_eq!(
            z.klein().unwrap(),
            Perplex::new(-1.0, 0.0),
            "-2 + h is in the left-sector!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::new(1.0, 2.0); // Up-Sector
        assert!(z.is_space_like(), "1 + 2h is space-like!");
        assert_ne!(z.arg(), 0.0, "1 + 2h has a non-zero argument!");
        assert_eq!(
            z.klein().unwrap(),
            Perplex::new(0.0, 1.0),
            "1 + 2h is in the up-sector!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::new(1.0, -2.0); // Down-Sector
        assert!(z.is_space_like(), "1 - 2h is space-like!");
        assert_ne!(z.arg(), 0.0, "1 - 2h has a non-zero argument!");
        assert_eq!(
            z.klein().unwrap(),
            Perplex::new(0.0, -1.0),
            "1 - 2h is in the down-sector!"
        );
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);

        let z = Perplex::cis(f64::PI() / 2.0);
        assert_abs_diff_eq!(z, Perplex::from(HyperbolicPolar::from(z)), epsilon = 0.0001);
    }

    fn polar_mul_test_loop(z: Perplex<f64>) {
        let polar = HyperbolicPolar::from(z);
        assert_abs_diff_eq!(Perplex::default(), Perplex::from(polar.pow(0)));
        assert_abs_diff_eq!(z, Perplex::from(polar.pow(1)), epsilon = 0.0001);
        assert_abs_diff_eq!(z * z, Perplex::from(polar.pow(2)), epsilon = 0.0001);
        assert_abs_diff_eq!(z * z * z, Perplex::from(polar.pow(3)), epsilon = 0.0001);
        assert_abs_diff_eq!(z * z * z * z, Perplex::from(polar.pow(4)), epsilon = 0.0001);
    }

    #[test]
    fn test_polar_multiplication() {
        let z = Perplex::new(1.0, 1.0); // Diagonal x=t
        polar_mul_test_loop(z);
        let z = Perplex::new(1.0, -1.0); // Diagonal x=-t
        polar_mul_test_loop(z);
        let z = Perplex::new(2.0, 1.0); // Right-Sector
        polar_mul_test_loop(z);
        polar_mul_test_loop(z.inv().unwrap());
        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        polar_mul_test_loop(z);
        polar_mul_test_loop(z.inv().unwrap());
        let z = Perplex::new(1.0, 2.0); // Up-Sector
        polar_mul_test_loop(z);
        polar_mul_test_loop(z.inv().unwrap());
        let z = Perplex::new(1.0, -2.0); // Down-Sector
        polar_mul_test_loop(z);
        polar_mul_test_loop(z.inv().unwrap());
        let z = Perplex::cis(f64::PI() / 2.0);
        polar_mul_test_loop(z);
        polar_mul_test_loop(z.inv().unwrap());
    }
    #[test]
    fn test_polar_sector() {
        let perplex = Perplex::new(1.0, 0.5);
        assert_eq!(perplex.sector(), HyperbolicSector::Right);
        let polar = perplex.polar();
        assert_eq!(polar.rho, f64::sqrt(0.75));
        assert_eq!(polar.theta, perplex.arg());
        assert_eq!(polar.sector, HyperbolicSector::Right);
    }
}
