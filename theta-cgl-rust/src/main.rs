#![allow(non_snake_case)]
use theta_cgl_rust::params::p254::Fp2;
use theta_cgl_rust::thp254::{Fq, ThetaPoint, CGL, ThetaPointDim2, CGL2};

fn dim1(X: Fp2, Z: Fp2) {
    println!("CGL dimension 1");
    let O0 = ThetaPoint::new(&X, &Z);

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

    let cgl = CGL::new(O0);
    let hash = cgl.hash(msg.to_vec());

    println!("{}", hash);
}

fn dim2(X: Fp2, Z: Fp2, U: Fp2, V: Fp2) {
    println!("CGL dimension 2");
    let O0 = ThetaPointDim2::new(&(X * U), &(X * V), &(Z * U), &(Z * V));

    let msg: [u8; 7] = [1, 1, 1, 0, 1, 1, 1];
    let cgl = CGL2::new(O0);
    let hash = cgl.hash(msg.to_vec(), 3);
    // let hash = cgl_hash_dim2(O0, msg.to_vec());

    println!("hash:");
    println!("{0}, {1}, {2}", &hash.0, &hash.1, &hash.2);
    println!("");
    println!("{:?}, {:?}, {:?}", &hash.0, &hash.1, &hash.2);
}


fn main() {
    // From sagemath
    // The corresponding theta null point for y^2 = x^3 + x
    //
    // a = i*12275542822304839524828636957574153606049258432202019246703196872283011715506 + 0
    // b = i + 17866357519039022340746304327512392032047517165206258904525681907470971174910
    //
    // Encoding function from Fp^2 element to bytes is available in the sage code
    let X_hex: &str = "0000000000000000000000000000000000000000000000000000000000000000b29164fbeafb402b03e1a2844c5f05e206d84f96f7287a2c76c92b6505b6231b";
    let Z_hex: &str = "feffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f270100000000000000000000000000000000000000000000000000000000000000";
    let (X, _) = Fq::decode(&hex::decode(X_hex).unwrap());
    let (Z, _) = Fq::decode(&hex::decode(Z_hex).unwrap());
    
    dim1(X, Z);
    // dim2(X, Z, X, Z);
}
