use criterion::{black_box, criterion_group, criterion_main, Criterion};
use num_traits::Pow;
use perplex_num::{HyperbolicPolar, Perplex, PerplexMatrixForm};
criterion_group!(benches, bench_multiplication);

criterion_main!(benches);
#[inline]
fn loop_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    let mut result = Perplex::new(1.0, 0.0);
    for _ in 0..exp {
        result *= z;
    }
    result
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
    group.bench_function("Perplex mul loop", |b| {
        b.iter(|| {
            let z = Perplex::new(0.123, 4.321);
            let exp = 50; // at 476 t becomes inf
            let _ = black_box(loop_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.bench_function("Matrix power", |b| {
        b.iter(|| {
            let z = Perplex::new(0.123, 4.321);
            let exp = 50; // at 476 t becomes inf
            let _ = black_box(matrix_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.bench_function("Polar form", |b| {
        b.iter(|| {
            let z = Perplex::new(0.123, 4.321);
            let exp = 50; // at 476 t becomes inf
            let _ = black_box(polar_multiplication(black_box(z), black_box(exp)));
        })
    });
    group.finish();
}
