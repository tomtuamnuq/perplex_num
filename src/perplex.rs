use std::ops::Neg;

use approx::AbsDiffEq;
use num_traits::float::FloatCore;
use num_traits::{Float, MulAdd, MulAddAssign, Num, NumAssign, One, Zero};

// The Mathematics of Minkowski Space-Time
// 4.1 Geometrical Representation of Hyperbolic Numbers

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
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

// TODO impl Display t + x h

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
    #[inline]
    fn default() -> Self {
        Self::new(T::one(), T::zero())
    }
}
impl<T: Copy + Num> Perplex<T> {
    /// Returns hyperbolic unit.
    #[inline]
    pub fn h() -> Self {
        Self::new(T::zero(), T::one())
    }
    #[inline]
    pub fn real(&self) -> T {
        self.t
    }
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
impl<T: Copy + Num + Neg<Output = T>> Perplex<T> {
    /// Returns the hyperbolic conjugate.
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.t, -self.x)
    }
    /// Returns the multiplicative inverse `1/self`, if it exists, or None if not.
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

    // TODO sqrt, cbrt, powf, log, powc, expf not existent for Perplex ?

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

impl<T: Copy + Num> From<T> for Perplex<T> {
    #[inline]
    fn from(t: T) -> Self {
        Self::new(t, T::zero())
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
