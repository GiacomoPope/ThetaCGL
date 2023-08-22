#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp254;

// sha256("Bristol 2023")
// static MSG: [u8; 256] = [
//     1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 
//     0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 
//     0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
//     1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 
//     0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 
//     0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 
//     0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 
//     1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1
// ];
static  MSG: [u8; 3] = [1, 1, 1];

fn dimension_one_254_example(X: thp254::Fq, Z: thp254::Fq) {
    println!("CGL dimension 1");
    let O0 = thp254::ThetaPoint::new(&X, &Z);
    let cgl = thp254::CGL::new(O0);
    let hash = cgl.hash(MSG.to_vec());

    println!("{}", hash);
}

fn dimension_one_127_example(X: thp127::Fq, Z: thp127::Fq) {
    println!("CGL dimension 1");
    let O0 = thp127::ThetaPoint::new(&X, &Z);
    let cgl = thp127::CGL::new(O0);
    let hash = cgl.hash(MSG.to_vec());

    println!("{}", hash);
}

// fn dimension_two_example(X: thp127::Fq, Z: thp127::Fq, U: thp127::Fq, V: thp127::Fq) {
//     println!("CGL dimension 2");
//     let O0 = thp127::ThetaPointDim2::new(&(X * U), &(X * V), &(Z * U), &(Z * V));

//     let msg: [u8; 7] = [1, 1, 1, 0, 1, 1, 1];
//     let cgl = thp127::CGL2::new(O0);
//     let hash: (thp127::Fq, thp127::Fq, thp127::Fq) = cgl.hash(msg.to_vec(), 3);
//     // let hash = cgl_hash_dim2(O0, msg.to_vec());

//     println!("hash:");
//     println!("{0}, {1}, {2}", &hash.0, &hash.1, &hash.2);
//     println!("");
//     println!("{:?}, {:?}, {:?}", &hash.0, &hash.1, &hash.2);
// }


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


    let X_hex_127: &str = "0000000000000000000000000000000000000000000000000100000000000000";
    let Z_hex_127: &str = "feffffffffffffffffffffffffffff7f01000000000000000000000000000000";
    let (X_127, _) = thp127::Fq::decode(&hex::decode(X_hex_127).unwrap());
    let (Z_127, _) = thp127::Fq::decode(&hex::decode(Z_hex_127).unwrap());

    dimension_one_127_example(X_127, Z_127);
}
