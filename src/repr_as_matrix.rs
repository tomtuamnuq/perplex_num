use super::perplex::Perplex;
use nalgebra::{Matrix2, RealField};
use num_traits::{One, Pow};
pub type PerplexMatrixForm<T> = Matrix2<T>; // TODO symmetric ? limit operations to the ones defined on a wrapping struct ?
impl<T: Copy + RealField> From<PerplexMatrixForm<T>> for Perplex<T> {
    fn from(m: PerplexMatrixForm<T>) -> Self {
        Self { t: m.m11, x: m.m12 }
    }
}
impl<T: Copy + RealField> From<Perplex<T>> for PerplexMatrixForm<T> {
    fn from(z: Perplex<T>) -> Self {
        Self::new(z.t, z.x, z.x, z.t)
    }
}

impl<T: Copy + RealField> Perplex<T> {
    #[inline]
    pub fn matrix_form(&self) -> PerplexMatrixForm<T> {
        let t = self.real();
        let x = self.hyperbolic();
        PerplexMatrixForm::new(t, x, x, t)
    }
    /// Raises ˋselfˋ to an unsigned integer power.
    #[inline]
    pub fn powu(&self, exp: u32) -> Self {
        Pow::pow(*self, exp)
    }

    /// Raises ˋselfˋ to a signed integer power.
    #[inline]
    pub fn powi(&self, exp: i32) -> Option<Self> {
        Pow::pow(*self, exp)
    }
}
// TODO impl for references ?
impl<T: Copy + RealField> Pow<u32> for Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn pow(self, exp: u32) -> Self::Output {
        // TODO benchmark if pow in matrix form is more performant
        match exp {
            0 => Perplex::one(),
            1 => self,
            2 => self * self,
            _ => {
                let m = PerplexMatrixForm::from(self);
                m.pow(exp).into()
            }
        }
    }
}

impl<T: Copy + RealField> Pow<i32> for Perplex<T> {
    type Output = Option<Perplex<T>>;
    #[inline]
    fn pow(self, exp: i32) -> Self::Output {
        if exp < 0 {
            PerplexMatrixForm::from(self)
                .try_inverse()
                .map(|m| m.pow(exp.wrapping_neg() as u32).into())
        } else {
            Some(Pow::pow(self, exp as u32))
        }
    }
}
