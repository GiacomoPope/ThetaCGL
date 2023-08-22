#![allow(non_snake_case)]
use theta_cgl_rust::params::p254::{Fp, Fp2};
use theta_cgl_rust::thp254::{Fq, ThetaPoint, cgl_hash};
use num_bigint::{BigInt, Sign};

fn pretty_print_Fp(v: &Fp) -> String{
    // Very stupid function, but it works...
    let v_bytes = v.encode();
    let v_big = BigInt::from_bytes_le(Sign::Plus, &v_bytes);
    v_big.to_string()
}

fn pretty_print_Fp2(v: &Fp2) -> String{
    // Very stupid function, but it works...
    let r = v.encode();

    let x0_bytes = &r[..Fp::ENCODED_LENGTH];
    let x1_bytes = &r[Fp::ENCODED_LENGTH..];

    let x0_big = BigInt::from_bytes_le(Sign::Plus, x0_bytes);
    let x1_big = BigInt::from_bytes_le(Sign::Plus, x1_bytes);

    format!("i*{} + {}", x1_big, x0_big)
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
    let hash = cgl_hash(O0, &msg);

    println!("{}", pretty_print_Fp2(&hash));


}
