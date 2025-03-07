#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_cgl_rust::thp5248;

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

fn two_radical_5248(c: &mut Criterion) {
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad2::new(block_size);
    c.bench_function(
        "CGL Hash (g=1) using 256-bit prime and two radical isogeny",
        |b| b.iter(|| cgl.hash(black_box(MSG.to_vec()))),
    );
}

fn four_radical_5248(c: &mut Criterion) {
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad4::new(block_size);
    c.bench_function(
        "CGL Hash (g=1) using 256-bit prime and four radical isogeny",
        |b| b.iter(|| cgl.hash(black_box(MSG.to_vec()))),
    );
}

fn eight_radical_5248(c: &mut Criterion) {
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad8::new(block_size);
    c.bench_function(
        "CGL Hash (g=1) using 256-bit prime and eight radical isogeny",
        |b| b.iter(|| cgl.hash(black_box(MSG.to_vec()))),
    );
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(15));
    targets = two_radical_5248, four_radical_5248, eight_radical_5248
}
criterion_main!(benches);
