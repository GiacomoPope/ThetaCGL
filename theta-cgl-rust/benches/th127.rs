#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_cgl_rust::thp127;

fn two_radical_127(c: &mut Criterion) {
    let X_hex_127: &str = "fdffffffffffffffffffffffffffff6b00000000000000000000000000000000";
    let Z_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let U_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let V_hex_127: &str = "00000000000000000000000000000000fdffffffffffffffffffffffffffff6b";

    let (X_127, _) = thp127::Fq::decode(&hex::decode(X_hex_127).unwrap());
    let (Z_127, _) = thp127::Fq::decode(&hex::decode(Z_hex_127).unwrap());
    let (U_127, _) = thp127::Fq::decode(&hex::decode(U_hex_127).unwrap());
    let (V_127, _) = thp127::Fq::decode(&hex::decode(V_hex_127).unwrap());

    let O0: thp127::ThetaPointDim2 = thp127::ThetaPointDim2::new(&X_127, &Z_127, &U_127, &V_127);
    let cgl = thp127::CGLDim2Rad2::new(O0);

    // sha256("Bristol 2023")
    let msg: [u8; 256] = [
        1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0,
        0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
        0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0,
        1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1,
        0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
    ];

    c.bench_function("CGL Hash: using p127", |b| {
        b.iter(|| cgl.hash(black_box(msg.to_vec()), black_box(3)))
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(10));
    targets = two_radical_127
}
criterion_main!(benches);
