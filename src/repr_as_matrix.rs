use super::perplex::Perplex;
use super::RealField;
use nalgebra::Matrix2;
pub type PerplexMatrixForm<T> = Matrix2<T>; // TODO symmetric ? limit operations to the ones defined on a wrapping struct ?
impl<T: RealField + Copy> From<PerplexMatrixForm<T>> for Perplex<T> {
    fn from(m: PerplexMatrixForm<T>) -> Self {
        Self { t: m.m11, x: m.m12 }
    }
}
impl<T: RealField + Copy> From<Perplex<T>> for PerplexMatrixForm<T> {
    fn from(z: Perplex<T>) -> Self {
        Self::new(z.t, z.x, z.x, z.t)
    }
}
