use std::ops::{Add, Div, Mul, Neg, Sub};

use num_traits::float::FloatCore;
use num_traits::{Float, Inv, MulAdd, Num, One, Zero};

// The Mathematics of Minkowski Space-Time
// 4.1 Geometrical Representation of Hyperbolic Numbers

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Perplex<T> {
    /// first coordinate t for time as real part
    pub t: T,
    /// second coordinate x for space as hyperbolic part
    pub x: T,
}

impl<T> Perplex<T> {
    /// Create a new Perplex number
    #[inline]
    pub fn new(t: T, x: T) -> Self {
        Self { t, x }
    }
}

// TODO impl Display t + x j

impl<T: Clone + Num> Perplex<T> {
    /// Returns hyperbolic unit.
    #[inline]
    pub fn h() -> Self {
        Self::new(T::zero(), T::one())
    }
    #[inline]
    pub fn real(&self) -> T {
        self.t.clone()
    }
    #[inline]
    pub fn hyperbolic(&self) -> T {
        self.x.clone()
    }
    /// Returns the squared distance D(z) in the hyperbolic plane.
    #[inline]
    pub fn squared_distance(&self) -> T {
        self.t.clone() * self.t.clone() - self.x.clone() * self.x.clone()
    }
    /// Multiplies ˋselfˋ by the scalar ˋfactorˋ.
    #[inline]
    pub fn scale(&self, factor: T) -> Self {
        Self::new(factor.clone() * self.t.clone(), factor * self.x.clone())
    }
}
impl<T: Clone + Num + PartialOrd> Perplex<T> {
    /// Checks if the perplex number is time-like, i.e., the squared distance is positive. See Properties of the Perplex Numbers in [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
    #[inline]
    pub fn is_time_like(&self) -> bool {
        self.squared_distance() > T::zero()
    }
    /// Checks if the perplex number is space-like, i.e., the squared distance is negative. See Properties of the Perplex Numbers in [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
    #[inline]
    pub fn is_space_like(&self) -> bool {
        self.squared_distance() < T::zero()
    }
    /// Checks if the perplex number is light-like, i.e., the squared distance is zero. See Properties of the Perplex Numbers in [Fundamental Theorems of Algebra for the Perplexes](https://doi.org/10.4169/074683409X475643)
    #[inline]
    pub fn is_light_like(&self) -> bool {
        self.squared_distance() == T::zero()
    }
}
impl<T: Clone + Num + Neg<Output = T>> Perplex<T> {
    /// Returns the hyperbolic conjugate.
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.t.clone(), -self.x.clone())
    }
    /// Returns the multiplicative inverse ˋ1/selfˋ, if it exists, or an error if not.
    #[inline]
    pub fn inv(&self) -> Self {
        let squared_distance = self.squared_distance();
        // TODO return Result - Err if zero
        Self::new(
            self.t.clone() / squared_distance.clone(),
            -self.x.clone() / squared_distance,
        )
    }
}

impl<T: Clone + Float> Perplex<T> {
    // TODO L1-Norm of a Perplex number given by Manhattan distance?

    /// Returns the modulus of ˋselfˋ.
    #[inline]
    pub fn modulus(self) -> T {
        let d_z = self.squared_distance();
        d_z.abs().sqrt()
    }
    /// Returns the norm (modulus) of ˋselfˋ.
    #[inline]
    pub fn norm(self) -> T {
        self.modulus()
    }
    /// Returns the magnitude (modulus) of ˋselfˋ.
    #[inline]
    pub fn magnitude(self) -> T {
        self.modulus()
    }

    // TODO polar coordinates: cis

    /// Calculate the argument of ˋselfˋ, if it exists. Formula is taken from Eq. 4.1.6 in Sec 4.1.1 Hyperbolic Exponential Function and Hyperbolic Polar Transformation in The Mathematics of Minkowski Space-Time.
    #[inline]
    pub fn arg(self) -> Option<T> {
        // TODO return Result
        let Self { t, x } = self;
        let (t_abs, x_abs) = (t.abs(), x.abs());
        if t_abs == x_abs {
            // light-like
            None
        } else if t_abs > x_abs {
            Some((x_abs / t_abs).atanh())
        } else {
            Some((t_abs / x_abs).atanh())
        }
    }
    /// Convert to hyperbolic polar form (rho, theta) such that ˋselfˋ = rho (cosh(theta) + h sinh(theta)), if theta (ˋself.arg()ˋ) exists.
    #[inline]
    pub fn to_polar(self) -> Option<(T, T)> {
        self.arg().map(|theta| (self.norm(), theta))
    }
    /// Convert from hyperbolic polar form (rho, theta) such that ˋSelfˋ = rho (cosh(theta) + h sinh(theta)). Formula is taken from Table 4.1 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6). The resulting ˋSelfˋ is in the right sector of the hyperbolic plane.
    #[inline]
    pub fn from_polar(rho: T, theta: T) -> Self {
        // TODO test if this provides the conversion given arg formula - how to place (x,y) into the 4 quadrants?
        Self::new(rho * theta.cosh(), rho * theta.sinh())
    }

    /// Computes the hyperbolic exponential function. Formula is taken from Sec 4.1.1 Hyperbolic Exponential Function and 7.4 The Elementary Functions of a Canonical Hyperbolic Variable in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn exp(self) -> Self {
        let Self { t, x } = self;
        let signum = if t.abs() > x.abs() {
            // self.is_time_like()
            t.signum()
        } else {
            x.signum()
        };
        let t_exp = signum * t.exp();
        Self::new(t_exp * x.cosh(), t_exp * x.sinh())
    }
    /// Computes the inverse of the hyperbolic exponential function, i.e., the natural logarithm. Formula is taken from Sec. 7.4 The Elementary Functions of a Canonical Hyperbolic Variable in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn ln(self) -> Option<Self> {
        let Self { t, x } = self;
        let squared_distance = t * t - x * x;
        let two = T::one() + T::one();
        if squared_distance == T::zero() {
            // light-like
            None
        } else if squared_distance > T::zero() {
            // self.is_time_like()
            let t_new = squared_distance.ln() / two;
            let x_new = (x / t).atanh();
            Some(Self::new(t_new, x_new))
        } else {
            // self.is_space_like()
            let t_new = squared_distance.neg().ln() / two;
            let x_new = (t / x).atanh();
            Some(Self::new(t_new, x_new))
        }
    }

    // TODO sqrt, cbrt, powf, log, powc, expf not existent for Perplex ?

    /// Computes the sinus (circular trigonometric) of ˋselfˋ. Formula is taken from Eq. 7.4.6 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn sin(self) -> Self {
        Self::new(self.t.sin() * self.x.cos(), self.t.cos() * self.x.sin())
    }
    /// Computes the cosinus (circular trigonometric) of ˋselfˋ. Formula is taken from Eq. 7.4.6 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn cos(self) -> Self {
        Self::new(self.t.cos() * self.x.cos(), self.t.sin() * self.x.sin())
    }
    /// Computes the sinh (hyperbolic trigonometric) of ˋselfˋ. Formula is taken from Eq. 7.4.5 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn sinh(self) -> Self {
        Self::new(self.t.sinh() * self.x.cosh(), self.t.cosh() * self.x.sinh())
    }
    /// Computes the cosh (hyperbolic trigonometric) of ˋselfˋ. Formula is taken from Eq. 7.4.5 in [The Mathematics of Minkowski Space-Time](https://doi.org/10.1007/978-3-7643-8614-6).
    #[inline]
    pub fn cosh(self) -> Self {
        Self::new(self.t.cosh() * self.x.cosh(), self.t.sinh() * self.x.sinh())
    }

    // TODO tan, asin, acos, atan, tanh, asinh, acosh, atanh ?
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

impl<T: Clone + Num> From<T> for Perplex<T> {
    #[inline]
    fn from(t: T) -> Self {
        Self::new(t, T::zero())
    }
}

impl<'a, T: Clone + Num> From<&'a T> for Perplex<T> {
    #[inline]
    fn from(t: &T) -> Self {
        From::from(t.clone())
    }
}

// Properties of the Perplex Numbers in Fundamental Theorems of Algebra for the Perplexes:
impl<T: Clone + Num> Add for Perplex<T> {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.t + rhs.t, self.x + rhs.x)
    }
}
impl<T: Clone + Num> Sub for Perplex<T> {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.t - rhs.t, self.x - rhs.x)
    }
}
impl<T: Clone + Num> Mul for Perplex<T> {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.t.clone() * rhs.t.clone() + self.x.clone() * rhs.x.clone(),
            rhs.t * self.x + self.t * rhs.x,
        )
    }
}
impl<T: Clone + Num + MulAdd<Output = T>> MulAdd<Perplex<T>> for Perplex<T> {
    type Output = Perplex<T>;

    #[inline]
    fn mul_add(self, other: Perplex<T>, add: Perplex<T>) -> Self {
        let t = self.t.clone() * other.t.clone() + self.x.clone() * other.x.clone() + add.t;
        let x = other.t * self.x + self.t * other.x + add.x;
        Self::new(t, x)
    }
}
impl<T: Clone + Num> Div for Perplex<T> {
    type Output = Option<Self>;
    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let Self { t: t2, x: x2 } = rhs;
        let norm_squared_2 = t2.clone() * t2.clone() - x2.clone() * x2.clone();
        if norm_squared_2 == T::zero() {
            // light-like
            None
        } else {
            let Self { t: t1, x: x1 } = self;
            let t_new =
                (t1.clone() * t2.clone() - x1.clone() * x2.clone()) / norm_squared_2.clone();
            let x_new = (t2 * x1 - t1 * x2) / norm_squared_2;
            Some(Self::new(t_new, x_new))
        }
    }
}

// TODO binary ops between references
// TODO assignment ops

impl<T: Clone + Num + Neg<Output = T>> Neg for Perplex<T> {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new(-self.t, -self.x)
    }
}

impl<'a, T: Clone + Num + Neg<Output = T>> Neg for &'a Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl<T: Clone + Num + Neg<Output = T>> Inv for Perplex<T> {
    type Output = Self;
    #[inline]
    fn inv(self) -> Self::Output {
        (&self).inv()
    }
}

impl<'a, T: Clone + Num + Neg<Output = T>> Inv for &'a Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn inv(self) -> Self::Output {
        self.inv()
    }
}

impl<T: Clone + Num> Add<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn add(self, other: T) -> Self::Output {
        Self::Output::new(self.t + other, self.x)
    }
}

impl<T: Clone + Num> Sub<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn sub(self, other: T) -> Self::Output {
        Self::Output::new(self.t - other, self.x)
    }
}

impl<T: Clone + Num> Mul<T> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn mul(self, other: T) -> Self::Output {
        Self::Output::new(self.t * other.clone(), self.x * other)
    }
}

impl<T: Clone + Num> Div<T> for Perplex<T> {
    type Output = Self;
    #[inline]
    fn div(self, other: T) -> Self::Output {
        Self::Output::new(self.t / other.clone(), self.x / other)
    }
}

// constants
impl<T: Clone + Num> Zero for Perplex<T> {
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

impl<T: Clone + Num> One for Perplex<T> {
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
