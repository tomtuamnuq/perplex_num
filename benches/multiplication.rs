use criterion::{black_box, criterion_group, criterion_main, Criterion};
use num_traits::Pow;
use perplex_num::{HyperbolicPolar, Perplex, PerplexMatrixForm};
criterion_group!(benches, bench_multiplication);
criterion_main!(benches);

const POW_EXP: u32 = 100;
const TIME: f64 = 0.123;
const SPACE: f64 = 4.321;
// for z=Perplex(TIME, SPACE) t becomes inf at z^476

#[inline]
fn loop_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    let mut result = Perplex::new(1.0, 0.0);
    for _ in 0..exp {
        result *= z;
    }
    result
}
#[inline]
fn squaring_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    z.powu(exp)
}
#[inline]
fn matrix_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    let m = PerplexMatrixForm::from(z);
    m.pow(exp).into()
}
#[inline]
fn polar_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    let polar = HyperbolicPolar::from(z);
    polar.pow(exp).into()
}
fn bench_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Multiplication");
    group.bench_function("Perplex mul naive loop", |b| {
        b.iter(|| {
            let z = Perplex::new(TIME, SPACE);
            let exp = POW_EXP;
            let _ = black_box(loop_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.bench_function("Perplex exponentiation by squaring", |b| {
        b.iter(|| {
            let z = Perplex::new(TIME, SPACE);
            let exp = POW_EXP;
            let _ = black_box(squaring_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.bench_function("Matrix power", |b| {
        b.iter(|| {
            let z = Perplex::new(TIME, SPACE);
            let exp = POW_EXP;
            let _ = black_box(matrix_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.bench_function("Polar form", |b| {
        b.iter(|| {
            let z = Perplex::new(TIME, SPACE);
            let exp = POW_EXP;
            let _ = black_box(polar_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.finish();
}
