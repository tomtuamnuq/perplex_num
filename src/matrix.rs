//! # Matrix Module
//!
//! This module is conditionally compiled only if the `matrix` feature is enabled. It provides an interface between the `Perplex` struct and the `nalgebra` crate for matrix operations.
//!
//! ## Overview
//! The `PerplexMatrixForm` type is a type alias for a 2x2 matrix from the `nalgebra` crate, representing the matrix form of a perplex number. This module includes conversions between the `Perplex` struct and its matrix representation.
//! The matrix representation of a perplex number is symmetric, with the real part on the diagonal and the hyperbolic part on the off-diagonal. This symmetry reflects the properties of perplex numbers regarding the hyperbolic unit.
//! Addition, multiplication, as well as inversion of perplex numbers correspond to the matrix operations.
//!
//! ## Example
//! ```rust
//! use perplex_num::Perplex;
//! let (z1, z2) = (Perplex::new(1.0, 2.0), Perplex::new(0.5, 0.1));
//! let (m1, m2) = (z1.as_matrix_form(), z2.as_matrix_form());
//! assert_eq!(z1 + z2, Perplex::from(m1 + m2), "Addition corresponds to matrix addition!");
//! assert_eq!(z1 * z2, Perplex::from(m1 * m2), "Multiplication corresponds to matrix multiplication!");
//! assert_eq!(z1.try_inverse().unwrap(), Perplex::from(m1.try_inverse().unwrap()), "Multiplicative inverse corresponds to matrix inverse!");
//! assert_eq!(z1.squared_distance(), m1.determinant(), "Squared distance corresponds to the determinant!");
//! ```

use super::perplex::Perplex;
use nalgebra::{Matrix2, RealField};

/// A type alias for a 2x2 matrix from `nalgebra`, representing a perplex number as a matrix.
pub type PerplexMatrixForm<T> = Matrix2<T>;

impl<T: Copy + RealField> From<PerplexMatrixForm<T>> for Perplex<T> {
    /// Converts a matrix form to a perplex number, assuming a symmetric matrix.
    fn from(m: PerplexMatrixForm<T>) -> Self {
        Self { t: m.m11, x: m.m12 }
    }
}

impl<T: Copy + RealField> From<Perplex<T>> for PerplexMatrixForm<T> {
    /// Returns the matrix form of the perplex number.
    fn from(z: Perplex<T>) -> Self {
        Self::new(z.t, z.x, z.x, z.t)
    }
}

impl<T: Copy + RealField> Perplex<T> {
    /// Creates a matrix form from a perplex number, resulting in a symmetric matrix.
    #[inline]
    pub fn as_matrix_form(&self) -> PerplexMatrixForm<T> {
        let t = self.real();
        let x = self.hyperbolic();
        PerplexMatrixForm::new(t, x, x, t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_matrix() {
        let (z1, z2) = (Perplex::new(1.0, 0.5), Perplex::new(-1.0, -2.0));
        let (m1, m2) = (z1.as_matrix_form(), PerplexMatrixForm::from(z2));
        assert_eq!(
            z1 + z2,
            Perplex::from(m1 + m2),
            "Matrix addition corresponds to addition of perplex numbers!"
        );
        assert_eq!(
            z1 * z2,
            Perplex::from(m1 * m2),
            "Matrix multiplication corresponds to multiplication of perplex numbers!"
        );
    }
}
