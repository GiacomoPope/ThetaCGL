#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp254;

// sha256("Bristol 2023")
static MSG: [u8; 256] = [
    1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 
    0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 
    0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 
    0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 
    0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 
    1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1
];

fn dimension_one_254_example(X: thp254::Fq, Z: thp254::Fq) {
    println!("CGL dimension 1");
    let O0 = thp254::ThetaPoint::new(&X, &Z);
    let cgl = thp254::CGL::new(O0);
    let hash = cgl.hash(MSG.to_vec());

    println!("{}", hash);
}

fn dimension_two_127_example(X: thp127::Fq, Z: thp127::Fq, U: thp127::Fq, V: thp127::Fq) {
    println!("CGL dimension 2");
    let O0: thp127::ThetaPointDim2 = thp127::ThetaPointDim2::new(&X, &Z, &U, &V);
    let cgl = thp127::CGL2::new(O0);
    let (h1, h2, h3) = cgl.hash(MSG.to_vec(), 3);

    println!("hash:");
    println!("{0}\n{1}\n{2}", h1, h2, h3);
}

fn main() {
    // From sagemath
    // The corresponding theta null point for y^2 = x^3 + x
    //
    // a = i*12275542822304839524828636957574153606049258432202019246703196872283011715506 + 0
    // b = i + 17866357519039022340746304327512392032047517165206258904525681907470971174910
    //
    // Encoding function from Fp^2 element to bytes is available in the sage code
    let X_hex_254: &str = "0000000000000000000000000000000000000000000000000000000000000000b29164fbeafb402b03e1a2844c5f05e206d84f96f7287a2c76c92b6505b6231b";
    let Z_hex_254: &str = "feffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f270100000000000000000000000000000000000000000000000000000000000000";
    let (X_254, _) = thp254::Fq::decode(&hex::decode(X_hex_254).unwrap());
    let (Z_254, _) = thp254::Fq::decode(&hex::decode(Z_hex_254).unwrap());
    
    dimension_one_254_example(X_254, Z_254);

    let X_hex_127: &str = "fdffffffffffffffffffffffffffff6b00000000000000000000000000000000";
    let Z_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let U_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let V_hex_127: &str = "00000000000000000000000000000000fdffffffffffffffffffffffffffff6b";

    let (X_127, _) = thp127::Fq::decode(&hex::decode(X_hex_127).unwrap());
    let (Z_127, _) = thp127::Fq::decode(&hex::decode(Z_hex_127).unwrap());
    let (U_127, _) = thp127::Fq::decode(&hex::decode(U_hex_127).unwrap());
    let (V_127, _) = thp127::Fq::decode(&hex::decode(V_hex_127).unwrap());

    dimension_two_127_example(X_127, Z_127, U_127, V_127);
}
