use super::Perplex;
use num_traits::{Float, One};

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub enum HyperbolicSector<T> {
    #[default]
    Right,
    Up,
    Left,
    Down,
    Diagonal(T),
}

impl<T: Copy + Float> From<Perplex<T>> for HyperbolicSector<T> {
    fn from(z: Perplex<T>) -> Self {
        let Perplex { t, x } = z;
        let (t_abs, x_abs) = (t.abs(), x.abs());
        if t_abs == x_abs {
            // light-like
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

/// Convert to hyperbolic polar form (rho, theta, klein) such that ˋselfˋ = klein rho (cosh(theta) + h sinh(theta)).
/// Convert from hyperbolic polar form (rho, theta, klein) such that ˋSelfˋ = rho (cosh(theta) + h sinh(theta)). Formula is taken from Tab. 1 and Appendix B in [Hyperbolic trigonometry in two-dimensional space-time geometry](https://doi.org/10.1393/ncb/i2003-10012-9). The resulting ˋSelfˋ is in the right sector of the hyperbolic plane.
// TODO describe light-like case
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct HyperbolicPolar<T> {
    pub rho: T,
    pub theta: T,
    pub sector: HyperbolicSector<T>,
}

impl<T: Copy + Float> From<Perplex<T>> for HyperbolicPolar<T> {
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
    /// Create a new Perplex with a given phase theta: `exp(h theta)`.
    #[inline]
    pub fn cis(theta: T) -> Self {
        Self::new(theta.cosh(), theta.sinh())
    }

    /// Calculate the hyperbolic argument of ˋselfˋ. Formula is taken from Eq. 4.1.6 in Sec 4.1.1 Hyperbolic Exponential Function and Hyperbolic Polar Transformation in The Mathematics of Minkowski Space-Time.
    /// Lines x=t and x=-t are mapped to ∞ and -∞, according to Sec. 4.1 in [New characterizations of the ring of the split-complex numbers and the field C of complex numbers and their comparative analyses](https://doi.org/10.48550/arXiv.2305.04586).
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

    /// Calculate the klein index of ˋselfˋ for space- or time-like numbers. Formula is taken from Tab. 1 and Appendix B in [Hyperbolic trigonometry in two-dimensional space-time geometry](https://doi.org/10.1393/ncb/i2003-10012-9).
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
}