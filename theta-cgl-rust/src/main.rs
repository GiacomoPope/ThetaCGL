#![allow(non_snake_case)]

use theta_cgl_rust::thp254::ThetaPoint;
use theta_cgl_rust::params::p254::{Fp, Fp2};
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

    format!("{} + i*{}", x0_big.to_string(), x1_big.to_string())
}


fn main() {
    println!("Hello, world!");
    
    let a = Fp::from_i64(1);
    let b = Fp::from_i64(2);
    let c = Fp2::new(&a, &b);



    println!("{}", pretty_print_Fp2(&d));

}
