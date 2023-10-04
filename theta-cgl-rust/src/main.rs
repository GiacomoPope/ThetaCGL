#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp254;
use theta_cgl_rust::thp921;

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

fn dimension_one_rad_2_254_example(X: thp254::Fq, Z: thp254::Fq) {
    println!("Computing using 2-radical isogenies...");
    let O0 = thp254::ThetaPointDim1::new(&X, &Z);
    let cgl = thp254::CGLDim1Rad2::new(O0);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected: &str = "i*15290852273938253774059106869425797954948837171016496465493101432860878286670 + 13712685150411413371193018930564063407253270277046650776399748023573137349444";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_254_example(X: thp254::Fq, Z: thp254::Fq) {
    println!("Computing using 4-radical isogenies...");
    let O0 = thp254::ThetaPointDim1::new(&X, &Z);
    let cgl = thp254::CGLDim1Rad4::new(O0, thp254::Fq::new(&thp254::Fp::from_i32(0), &thp254::Fp::from_i32(1)));
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*17632232527346586858516803369370617443794677771378157522072067056904912338941 + 10075930966783318325990846713272328106171814976421617495054657906335886232969";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_2_921_example(X: thp921::Fq, Z: thp921::Fq) {
    println!("Computing using 2-radical isogenies...");
    let O0 = thp921::ThetaPointDim1::new(&X, &Z);
    let cgl = thp921::CGLDim1Rad2::new(O0);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected: &str = "i*39230057073311561979904713106661291660670794873352665774512435740425423658735 + 28281974817376390026165203588145415349405195023788259912005133424119065028093";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_921_example(X: thp921::Fq, Z: thp921::Fq) {
    println!("Computing using 4-radical isogenies...");
    let O0 = thp921::ThetaPointDim1::new(&X, &Z);
    let cgl = thp921::CGLDim1Rad4::new(O0, thp921::Fq::new(&thp921::Fp::from_i32(0), &thp921::Fp::from_i32(1)));
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*44254100953703562411049248400500098463853458943478030726661718620302339064539 + 539165798661926336559847596085035641527652935389270773552568291982882834598";
    println!("SageMath: {}", expected);
}

fn dimension_two_rad_2_127_example(X: thp127::Fq, Z: thp127::Fq, U: thp127::Fq, V: thp127::Fq) {
    println!("Computing using 2-radical isogenies...");
    let O0: thp127::ThetaPointDim2 = thp127::ThetaPointDim2::new(&X, &Z, &U, &V);
    let cgl = thp127::CGLDim2Rad2::new(O0);
    let (h1, h2, h3) = cgl.hash(MSG.to_vec(), 3);

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let exp1 = "i*91066962968548325954088205255508107847 + 120247548916298594825558313821946927041";
    let exp2 = "i*108091552431835511670716536575836453590 + 126978242929844491058319474121908139051";
    let exp3 = "i*93365410772160919332300248626962422489 + 25264642909927975483848712679055430268";

    println!("Sage hash: {}", exp1);
    println!("           {0}\n           {1}", exp2, exp3);
}

fn main() {    
    println!("================================================================================");
    println!("                  Dimension One CGL with p = 79*2^247 - 1"                       );
    println!("================================================================================");

    let X_hex_254: &str = "0000000000000000000000000000000000000000000000000000000000000000b29164fbeafb402b03e1a2844c5f05e206d84f96f7287a2c76c92b6505b6231b";
    let Z_hex_254: &str = "feffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f270100000000000000000000000000000000000000000000000000000000000000";
    let (X_254, _) = thp254::Fq::decode(&hex::decode(X_hex_254).unwrap());
    let (Z_254, _) = thp254::Fq::decode(&hex::decode(Z_hex_254).unwrap());

    dimension_one_rad_2_254_example(X_254, Z_254);
    dimension_one_rad_4_254_example(X_254, Z_254);

    
    println!("================================================================================");
    println!("                  Dimension One CGL with p = 2^255 - 921"                        );
    println!("================================================================================");

    let X_hex_921: &str = "0000000000000000000000000000000000000000000000000000000000000000c4c0c9cff0aa2eff0b84b935338c4131d4260a06e1fe83fcfb118559864c7a50";
    let Z_hex_921: &str = "66fcffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f0100000000000000000000000000000000000000000000000000000000000000";
    let (X_921, _) = thp921::Fq::decode(&hex::decode(X_hex_921).unwrap());
    let (Z_921, _) = thp921::Fq::decode(&hex::decode(Z_hex_921).unwrap());

    dimension_one_rad_2_921_example(X_921, Z_921);
    dimension_one_rad_4_921_example(X_921, Z_921);

    println!("================================================================================");
    println!("                  Dimension One CGL with p = 27*2^122 - 1"                       );
    println!("================================================================================");

    let X_hex_127: &str = "fdffffffffffffffffffffffffffff6b00000000000000000000000000000000";
    let Z_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let U_hex_127: &str = "1fc93e85eba36a2d4d49a011ce720f421fc93e85eba36a2d4d49a011ce720f42";
    let V_hex_127: &str = "00000000000000000000000000000000fdffffffffffffffffffffffffffff6b";

    let (X_127, _) = thp127::Fq::decode(&hex::decode(X_hex_127).unwrap());
    let (Z_127, _) = thp127::Fq::decode(&hex::decode(Z_hex_127).unwrap());
    let (U_127, _) = thp127::Fq::decode(&hex::decode(U_hex_127).unwrap());
    let (V_127, _) = thp127::Fq::decode(&hex::decode(V_hex_127).unwrap());
    
    dimension_two_rad_2_127_example(X_127, Z_127, U_127, V_127);
}
