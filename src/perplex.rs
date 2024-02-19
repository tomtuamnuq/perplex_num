use super::repr_as_matrix::PerplexMatrixForm;
use super::RealField;
// The Mathematics of Minkowski Space-Time TODO add DOI
// 4.1 Geometrical Representation of Hyperbolic Numbers
// first coordinate t for time as real part, second x for space as hyperbolic part
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct Perplex<T> {
    pub t: T,
    pub x: T,
}
// TODO impl Display t + x j
impl<T: Copy> Perplex<T> {
    pub fn new(t: T, x: T) -> Self {
        Self { t, x }
    }
    pub fn real(&self) -> T {
        self.t
    }
    pub fn hyperbolic(&self) -> T {
        self.x
    }
    pub fn representation(&self) -> PerplexMatrixForm<T> {
        let t = self.real();
        let x = self.hyperbolic();
        PerplexMatrixForm::new(t, x, x, t)
    }
}

impl<T: RealField + Copy> Perplex<T> {
    // Properties of the Perplex Numbers in Fundamental Theorems of Algebra for the Perplexes TODO add DOI
    // also called split-complex numbers or spacetime numbers
    // also: Sec 4.1 Geometrical Representations of Hyperbolic Numbers in The Mathematics of Minkowski Space-Time
    pub fn is_time_like(&self) -> bool {
        self.real().abs() > self.hyperbolic().abs()
    }
    pub fn is_space_like(&self) -> bool {
        self.real().abs() < self.hyperbolic().abs()
    }
    pub fn is_light_like(&self) -> bool {
        self.real().abs() == self.hyperbolic().abs()
    }
    pub fn magnitude(&self) -> Option<T> {
        if self.is_space_like() {
            None
        } else {
            let t = self.real();
            let x = self.hyperbolic();
            let norm = (t.powi(2) - x.powi(2)).sqrt();
            Some(norm)
        }
    }
    pub fn conjugate(self) -> Self {
        Self::new(self.t, -self.x)
    }
    pub fn inv(self) -> Option<Self> {
        let Self { t, x } = self;
        let norm_squared = t.powi(2) - x.powi(2);
        if norm_squared == T::zero() {
            // self.is_light_like()
            None
        } else {
            // 1/z = z.conjugate() / (z * z.conjugate)
            Some(Self::new(t / norm_squared, -x / norm_squared))
        }
    }

    // First Table in Hybrid Numbers Mustafa Özdemir TODO add DOI
    // also: Sec 4.1.1 Hyperbolic Exponential Function and Hyperbolic Polar Transformation in The Mathematics of Minkowski Space-Time
    pub fn argument(&self) -> Option<T> {
        let Self { t, x } = *self;
        let norm_squared = t.powi(2) - x.powi(2);
        if norm_squared == T::zero() {
            // light-like
            None
        } else {
            let norm = norm_squared.sqrt();
            let t_x = (t + x).abs();
            Some((t_x / norm).ln())
        }
    }

    // Sec 4.1.1 Hyperbolic Exponential Function and 7.4 The Elementary Functions of a Canonical Hyperbolic Variable in The Mathematics of Minkowski Space-Time
    pub fn exp(self) -> Self {
        let signum = if self.is_space_like() {
            self.x.signum()
        } else {
            self.t.signum()
        };
        let t_exp = signum * self.t.exp();
        Self::new(t_exp * self.x.cosh(), t_exp * self.x.sinh())
    }
    pub fn ln(self) -> Option<Self> {
        let Self { t, x } = self;
        let norm_squared = t.powi(2) - x.powi(2);
        let two: T = T::from_u8(2).unwrap();
        if norm_squared > T::zero() {
            // self.is_time_like()
            let t_new = norm_squared.ln() / two;
            let x_new = (x / t).atanh();
            Some(Self::new(t_new, x_new))
        } else if norm_squared < T::zero() {
            // self.is_space_like()
            let t_new = norm_squared.neg().ln() / two;
            let x_new = (t / x).atanh();
            Some(Self::new(t_new, x_new))
        } else {
            None
        }
    }
    pub fn cosh(self) -> Self {
        Self::new(self.t.cosh() * self.x.cosh(), self.t.sinh() * self.x.sinh())
    }
    pub fn sinh(self) -> Self {
        Self::new(self.t.sinh() * self.x.cosh(), self.t.cosh() * self.x.sinh())
    }
    pub fn cos(self) -> Self {
        Self::new(self.t.cos() * self.x.cos(), self.t.sin() * self.x.sin())
    }
    pub fn sin(self) -> Self {
        Self::new(self.t.sin() * self.x.cos(), self.t.cos() * self.x.sin())
    }

    // TODO impl traits from nalgebra?
    // TODO abs() of perplex number just elementwise?
    pub fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }
    pub fn one() -> Self {
        Self::new(T::one(), T::zero())
    }
    pub fn pow(self, n: u32) -> Self {
        // TODO benchmark if pow in matrix form is more performant
        match n {
            0 => Self::one(),
            1 => self,
            2 => self * self,
            _ => {
                let m = PerplexMatrixForm::from(self);
                m.pow(n).into()
            }
        }
    }
}

// Properties of the Perplex Numbers in Fundamental Theorems of Algebra for the Perplexes TODO add DOI
impl<T: RealField + Copy> std::ops::Add for Perplex<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.t + rhs.t, self.x + rhs.x)
    }
}
impl<T: RealField + Copy> std::ops::Sub for Perplex<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.t - rhs.t, self.x - rhs.x)
    }
}
impl<T: RealField + Copy> std::ops::Mul for Perplex<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.t * rhs.t + self.x * rhs.x,
            rhs.t * self.x + self.t * rhs.x,
        )
    }
}
impl<T: RealField + Copy> std::ops::Div for Perplex<T> {
    type Output = Option<Self>;
    fn div(self, rhs: Self) -> Self::Output {
        let Self { t: t2, x: x2 } = rhs;
        let norm_squared_2 = t2.powi(2) - x2.powi(2);
        if norm_squared_2 == T::zero() {
            None
        } else {
            let Self { t: t1, x: x1 } = self;
            let t_new = (t1 * t2 - x1 * x2) / norm_squared_2;
            let x_new = (t2 * x1 - t1 * x2) / norm_squared_2;
            Some(Self::new(t_new, x_new))
        }
    }
}
impl<T: RealField + Copy> std::ops::Mul<T> for Perplex<T> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        Self::new(rhs * self.t, rhs * self.x)
    }
}
