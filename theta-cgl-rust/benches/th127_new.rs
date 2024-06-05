#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_cgl_rust::thp127new;

// sha256("Bristol 2023")
const MSG: [u8; 256] = [
    1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0,
    0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,
    0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1,
    0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
];

fn two_radical_127(c: &mut Criterion) {
    let cgl = thp127new::CGLDim2Rad2::new();

    c.bench_function(
        "CGL Hash: using p127 and two radical isogeny and new arithmetic",
        |b| b.iter(|| cgl.hash(black_box(MSG.to_vec()))),
    );
}

fn four_radical_127(c: &mut Criterion) {
    let cgl = thp127new::CGLDim2Rad4::new();

    c.bench_function(
        "CGL Hash: using p127 and four radical isogeny and new arithmetic",
        |b| b.iter(|| cgl.hash(black_box(MSG.to_vec()))),
    );
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(15));
    targets = two_radical_127, four_radical_127
}
criterion_main!(benches);
