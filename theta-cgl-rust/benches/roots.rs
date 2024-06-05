mod fpcore;
use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_cgl_rust::fields::Fp127Old::Fp;

fn sqrt(c: &mut Criterion) {
    let a = Fp::new([0x0123456789abcdef, 0xfedcba9876543210]);
    let aa = a.square();

    c.bench_function("Square-root", |b| b.iter(|| aa.sqrt()));
}

fn fourth_root(c: &mut Criterion) {
    let a = Fp::new([0x0123456789abcdef, 0xfedcba9876543210]);
    let aa = a.square().square();

    c.bench_function("Fourth-root", |b| b.iter(|| aa.sqrt()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(10));
    targets = sqrt, fourth_root
}
criterion_main!(benches);
