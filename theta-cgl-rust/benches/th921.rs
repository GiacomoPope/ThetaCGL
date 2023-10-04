#![allow(non_snake_case)]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use theta_cgl_rust::thp921;



fn criterion_benchmark(c: &mut Criterion) {
    let X_hex: &str = "0000000000000000000000000000000000000000000000000000000000000000c4c0c9cff0aa2eff0b84b935338c4131d4260a06e1fe83fcfb118559864c7a50";
    let Z_hex: &str = "66fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f0100000000000000000000000000000000000000000000000000000000000000";
    let (X, _) = thp921::Fq::decode(&hex::decode(X_hex).unwrap());
    let (Z, _) = thp921::Fq::decode(&hex::decode(Z_hex).unwrap());
    let O0 = thp921::ThetaPointDim1::new(&X, &Z);
    let cgl = thp921::CGLDim1Rad2::new(O0);
    
    // sha256("Bristol 2023")
    let msg: [u8; 256] = [
        1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 
        0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 
        0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
        1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 
        0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 
        0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 
        0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 
        1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1
    ];

    c.bench_function(
        "CGL Hash: using p921", 
        |b| b.iter(|| cgl.hash(
             black_box(msg.to_vec())
        )
    ));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
