use super::perplex::Perplex;
use nalgebra::{Matrix2, RealField};
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
    pub fn as_matrix_form(&self) -> PerplexMatrixForm<T> {
        let t = self.real();
        let x = self.hyperbolic();
        PerplexMatrixForm::new(t, x, x, t)
    }
}
