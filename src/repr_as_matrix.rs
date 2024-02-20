use super::perplex::Perplex;
use nalgebra::{Matrix2, RealField};
use num_traits::{Num, One, Pow};
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
        Pow::pow(self, exp)
    }

    /// Raises ˋselfˋ to a signed integer power.
    #[inline]
    pub fn powi(&self, exp: i32) -> Option<Self> {
        Pow::pow(self, exp)
    }
}
impl<'a, T: Copy + RealField> Pow<u32> for &'a Perplex<T> {
    type Output = Perplex<T>;
    #[inline]
    fn pow(self, exp: u32) -> Self::Output {
        // TODO benchmark if pow in matrix form is more performant
        match exp {
            0 => Perplex::one(),
            1 => self.clone(),
            2 => self.clone() * self.clone(), // TODO impl for references
            _ => {
                let m = PerplexMatrixForm::from(self.clone());
                m.pow(exp).into()
            }
        }
    }
}

impl<'a, 'b, T: Copy + RealField> Pow<&'b u32> for &'a Perplex<T> {
    type Output = Perplex<T>;

    #[inline]
    fn pow(self, exp: &u32) -> Self::Output {
        self.pow(*exp)
    }
}

impl<'a, T: Copy + RealField> Pow<i32> for &'a Perplex<T> {
    type Output = Option<Perplex<T>>;
    #[inline]
    fn pow(self, exp: i32) -> Self::Output {
        if exp < 0 {
            PerplexMatrixForm::from(self.clone())
                .try_inverse()
                .map(|m| m.pow(exp.wrapping_neg() as u32).into())
        } else {
            Some(Pow::pow(self, exp as u32))
        }
    }
}

impl<'a, 'b, T: Copy + RealField> Pow<&'b i32> for &'a Perplex<T> {
    type Output = Option<Perplex<T>>;

    #[inline]
    fn pow(self, exp: &i32) -> Self::Output {
        self.pow(*exp)
    }
}
