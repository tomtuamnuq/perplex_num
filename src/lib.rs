mod binary_ops;
mod perplex;
mod polar;
mod repr_as_matrix;
pub use perplex::Perplex;
pub use polar::{HyperbolicPolar, HyperbolicSector};
pub use repr_as_matrix::PerplexMatrixForm;
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use num_traits::*;

    #[test]
    fn test_norm() {
        let z = Perplex::new(2.0, -1.0);
        assert!(z.is_time_like());
        assert_eq!(z.modulus(), f64::sqrt(3.0), "2 - h has a norm of √3");
        let z = Perplex::new(1.0, -1.0);
        assert!(z.is_light_like());
        assert_eq!(z.magnitude(), f64::zero(), "1 - h has a norm of zero");
        let z = Perplex::new(-1.0, 2.0);
        assert!(z.is_space_like());
        assert_eq!(z.norm(), f64::sqrt(3.0), "-1 + 2h has a norm of √3");
        assert_eq!(z.l1_norm(), 3.0, "-1 + 2h has a l1 norm of 3");
        assert_eq!(z.l2_norm(), f64::sqrt(5.0), "-1 + 2h has a l2 norm of √5");
        assert_eq!(z.max_norm(), 2.0, "-1 + 2h has a max norm of 2");
    }
    #[test]
    fn test_inv() {
        let z = Perplex::new(2.0, -1.0);
        let inv_result = z.inv();
        assert!(inv_result.is_some(), "2 - h is invertible!");
        assert_eq!(inv_result.unwrap(), Perplex::new(2.0 / 3.0, 1.0 / 3.0));
    }

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
    fn test_logarithm_exponential() {
        let z = Perplex::new(2.0, 1.0); // Right-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for time-like 2 + h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for time-like -2 + h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(1.0, 2.0); // Up-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for space-like 1 + 2h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);

        let z = Perplex::new(1.0, -2.0); // Down-Sector
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "Natural logarithm is defined for space-like 1 - 2h!"
        );
        let z_ln_exp = ln_result.unwrap().exp();
        assert_abs_diff_eq!(z_ln_exp, z);
    }
    #[test]
    fn test_exponential_logarithm() {
        let z = Perplex::new(2.0, 1.0); // Right-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(-2.0, 1.0); // Left-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(1.0, 2.0); // Up-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
        let z = Perplex::new(1.0, -2.0); // Down-Sector
        assert_abs_diff_eq!(z.exp().ln().unwrap(), z, epsilon = 0.00001);
    }

    #[test]
    fn test_trigonometric() {
        let pi = f64::PI();
        let z = Perplex::new(pi, pi / 2.0).sin();
        assert_abs_diff_eq!(z, Perplex::new(0.0, -1.0));
        let zero: Perplex<f64> = Perplex::zero();
        assert_abs_diff_eq!(zero.sinh(), zero);
    }

    #[test]
    fn test_add() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 + one + zero,
            Perplex::new(2.0, 2.0),
            "Componentwise addition!"
        );
        assert_eq!(
            z1 + z1.conj(),
            Perplex::new(2.0, 0.0),
            "Addition of conjugate zeros the hyperbolic part!"
        );
        let mut z2 = Perplex::new(-3.0, 2.0);
        let z12 = z1 + z2;
        z2 += z1;
        assert_eq!(z12, z2, "AddAssign yields same result as Add!");
    }
    #[test]
    fn test_sub() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 - one - zero,
            Perplex::new(0.0, 2.0),
            "Componentwise subtraction!"
        );
        assert_eq!(
            z1 - z1.conj(),
            Perplex::new(0.0, 4.0),
            "Subtraction of conjugate doubles the hyperbolic part!"
        );
        let mut z2 = Perplex::new(-3.0, 2.0);
        let z12 = z2 - z1;
        z2 -= z1;
        assert_eq!(z12, z2, "SubAssign yields same result as Sub!");
    }
    #[test]
    fn test_mul() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 * one,
            z1,
            "Neutral element of multiplication yields same element!"
        );
        assert_eq!(z1 * zero, zero, "Neutral element of addition yields zero!");
        let mut z2 = Perplex::new(-1.0, 2.0);
        let z12 = z1 * z2;
        z2 *= z1;
        assert_eq!(z2, Perplex::new(3.0, 0.0), "Multiplication formula!");
        assert_eq!(z12, z2, "MulAssign yields same result as Mul!");
    }
    #[test]
    fn test_div() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            (z1 / one).unwrap(),
            z1,
            "Division of neutral element of multiplication yields same element!"
        );
        assert!(
            (z1 / zero).is_none(),
            "Division of neutral element of addition yields none!"
        );
        let z2 = Perplex::new(-1.0, 2.0);
        let mut z12 = z1 * z2;
        let div_result = z12 / z2;
        assert!(
            div_result.is_some(),
            "Division of product by multiplier is valid!"
        );
        assert_eq!(
            div_result.unwrap(),
            z1,
            "Division of product by multiplier gives multiplicand."
        );
        z12 /= z2;
        assert_eq!(z12, z1, "DivAssign yields same result as Div!");

        let z2 = Perplex::new(-1.0, 1.0);
        let z12 = z1 * z2;
        assert_eq!(z12, Perplex::new(1.0, -1.0), "Multiplication formula!");
        assert!(z2.is_light_like(), "-1 + j is light-like!");
        assert!(
            (z12 / z2).is_none(),
            "Division is not defined for light-like numbers!"
        );
    }
    #[test]
    fn test_scalar() {
        let z1 = Perplex::new(1.0, 2.0);
        assert_eq!(
            z1 + 2.0,
            Perplex::new(3.0, 2.0),
            "Addition of scalar only on time component!"
        );
        assert_eq!(
            z1 - 2.0,
            Perplex::new(-1.0, 2.0),
            "Subtraction of scalar only on time component!"
        );
        assert_eq!(
            z1 * 2.0,
            Perplex::new(2.0, 4.0),
            "Componentwise scalar multiplication!"
        );
        assert_eq!(
            z1 / 2.0,
            Perplex::new(0.5, 1.0),
            "Componentwise scalar division!"
        );
    }

    #[test]
    fn test_power_u32() {
        let z = Perplex::new(1.0, -1.0);
        assert_eq!(
            z * z,
            Perplex::new(2.0, -2.0),
            "Multiplication with itself!"
        );
        assert_eq!(
            z.powu(2),
            Perplex::new(2.0, -2.0),
            "Power 2 yields multiplication with itself!"
        );
        assert_eq!(
            z.powu(3),
            Perplex::new(4.0, -4.0),
            "Power 3 multiplication result!"
        );
    }
    #[test]
    fn test_power_i32() {
        let z = Perplex::new(1.0, -1.0);
        assert!(z.powi(-2).is_none(), " 1 - h is not invertibe!");
        let z = Perplex::new(2.0, 1.0);
        let z_inv = z.try_inverse().unwrap();
        assert_eq!(
            z_inv.powi(2).unwrap(), // (2/3 - h / 3)^2
            Perplex::new(5.0 / 9.0, -4.0 / 9.0),
            "Multiplication of inverse with itself!"
        );
        assert_eq!(
            z.powi(-2).unwrap(),
            Perplex::new(5.0 / 9.0, -4.0 / 9.0),
            "Power -2 yields multiplication of inverse with itself!"
        );
        assert_eq!(
            z.powi(-3).unwrap(),
            Perplex::new(14.0 / 27.0, -13.0 / 27.0),
            "Power -3 multiplication result!"
        );
    }
}
