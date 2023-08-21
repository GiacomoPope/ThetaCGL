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

    let x0_big = BigInt::from_bytes_le(Sign::Plus, &x0_bytes);
    let x1_big = BigInt::from_bytes_le(Sign::Plus, &x1_bytes);

    format!("i*{} + {}", x1_big.to_string(), x0_big.to_string())
}


fn main() {
    // TODO: get proper null point
    let a = Fp::from_i64(1);
    let b = Fp::from_i64(2);
    let X: Fq = Fq::new(&a, &b);
    let Z: Fq = Fq::new(&b, &a);

    let O0 = ThetaPoint{
        X, Z
    };

    // let msg: [u8; 7] = [1, 1, 1, 0, 1, 1, 1];
    let msg: [u8; 1] = [1];

    let hash = cgl_hash(O0, &msg);

    println!("{}", pretty_print_Fp2(&hash));


}
