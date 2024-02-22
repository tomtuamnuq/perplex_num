use std::ops::MulAssign;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use perplex_num::Perplex;

criterion_group!(benches, bench_multiplication);

criterion_main!(benches);
#[inline]
fn loop_multiplication(z: Perplex<f64>, exp: u32) -> Perplex<f64> {
    match exp {
        0 => Perplex::new(1.0, 0.0),
        1 => z,
        _ => {
            let mut z_pow = black_box(z);
            for _ in 1..exp {
                black_box(MulAssign::mul_assign(black_box(&mut z_pow), black_box(z)));
            }
            black_box(z_pow)
        }
    }
}
fn bench_multiplication(c: &mut Criterion) {
    let mut group = c.benchmark_group("Multiplication");
    group.bench_function("Perplex mul loop", |b| {
        b.iter(|| {
            black_box(|| {
                let z = black_box(Perplex::new(0.123, 4.321));
                let exp = black_box(1000);
                black_box(loop_multiplication(black_box(z), black_box(exp)));
            })
        })
    });
    group.bench_function("Matrix power", |b| {
        b.iter(|| {
            let z = black_box(Perplex::new(0.123, 4.321));
            let exp = black_box(1000);
            black_box(black_box(z).powu(black_box(exp)));
        })
    });
    group.finish();
}
