mod perplex;
mod repr_as_matrix;
pub use perplex::Perplex;
pub use repr_as_matrix::PerplexMatrixForm;
#[cfg(test)]
mod tests {
    use super::*;
    fn approx_equal(z1: Perplex<f64>, z2: Perplex<f64>) -> bool {
        (z1.t - z2.t).abs() < 0.00001 && (z1.x - z2.x).abs() < 0.00001
    } // TODO impl AbsDiffEq
    #[test]
    fn test_std_ops() {
        let z1 = Perplex::new(1.0, 2.0);
        let one = Perplex::one();
        let zero = Perplex::zero();
        assert_eq!(
            z1 + one + zero,
            Perplex::new(2.0, 2.0),
            "Componentwise addition!"
        );
        assert_eq!(
            z1 + z1.conjugate(),
            Perplex::new(2.0, 0.0),
            "Addition of conjugate zeros the hyperbolic part!"
        );
        let h = Perplex::new(0.0, 1.0);
        assert_eq!(z1 - one - h, h, "Componentwise subtraction!");
        assert_eq!(
            z1 * one,
            z1,
            "Neutral element of multiplication yields same element!"
        );
        assert_eq!(z1 * zero, zero, "Neutral element of addition yields zero!");
        let z2 = Perplex::new(-1.0, 2.0);
        let z12 = z1 * z2;
        assert_eq!(z12, Perplex::new(3.0, 0.0), "Multiplication formula!");
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
        let z2 = Perplex::new(-1.0, 1.0);
        let z12 = z1 * z2;
        assert_eq!(z12, Perplex::new(1.0, -1.0), "Multiplication formula!");
        assert!(z2.is_light_like(), "-1 + j is light-like!");
        assert!(
            (z12 / z2).is_none(),
            "Division is not defined for light-like numbers!"
        );
        assert_eq!(
            z1 * 2.0,
            Perplex::new(2.0, 4.0),
            "Componentwise scalar multiplication!"
        );
    }
    #[test]
    fn test_trigonometric() {
        let pi = f64::pi();
        let z = Perplex::new(pi, pi / 2.0).sin();
        assert!(approx_equal(z, Perplex::new(0.0, -1.0)), "Sinus Formula!");
        let zero: Perplex<f64> = Perplex::zero();
        assert!(
            approx_equal(zero.sinh(), zero),
            "Sinus Hyperbolicus is zero at 0.0"
        );
        let z = Perplex::new(1.0, 2.0);
        dbg!(z.exp());
        assert!(
            approx_equal(z.exp().ln().unwrap(), z),
            "Logarithm inverts exponential function!"
        );
        assert!(z.is_space_like(), "1 + 2j is space like!");
        let ln_result = z.ln();
        assert!(
            ln_result.is_some(),
            "1 + 2j is space like and therefore is ln defined!"
        );
        dbg!(ln_result.unwrap());
        let z_ln_exp = ln_result.unwrap().exp();
        dbg!(z_ln_exp);
        assert!(
            approx_equal(z_ln_exp, z),
            "Exponential inverts logarithmic function!"
        );
    }
    #[test]
    fn test_norm() {
        let z = Perplex::new(2.0, -1.0);
        assert!(z.is_time_like());
        let inv_result = z.inv();
        assert!(inv_result.is_some(), "2 - j is invertible!");
        assert_eq!(z.magnitude(), Some(f64::sqrt(3.0)));
        assert_eq!(inv_result.unwrap(), Perplex::new(2.0 / 3.0, 1.0 / 3.0));
    }
    #[test]
    fn test_argument() {
        let z = Perplex::new(1.0, -1.0);
        assert!(z.is_light_like(), "1 - j is light-like!");
        assert!(
            z.argument().is_none(),
            "Argument is not defined for light-like numbers!"
        );
        let z = Perplex::new(1.0, 2.0);
        assert!(z.is_space_like(), "2 - j is space-like!");
        assert_ne!(z.argument().unwrap(), 0.0, "2 - j has a non-zero argument!");
    }
    #[test]
    fn test_power() {
        let z = Perplex::new(1.0, -1.0);
        assert_eq!(
            z * z,
            Perplex::new(2.0, -2.0),
            "Multiplication with itself!"
        );
        assert_eq!(
            z.pow(2),
            Perplex::new(2.0, -2.0),
            "Power 2 yields multiplication with itself!"
        );
        assert_eq!(
            z.pow(3),
            Perplex::new(4.0, -4.0),
            "Power 3 multiplication result is correct!"
        );
    }
}
