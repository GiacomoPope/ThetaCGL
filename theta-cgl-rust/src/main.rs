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
    1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
];

fn dimension_one_rad_2_254_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp254::CGLDim1Rad2::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected: &str = "i*15290852273938253774059106869425797954948837171016496465493101432860878286670 + 13712685150411413371193018930564063407253270277046650776399748023573137349444";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_254_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl = thp254::CGLDim1Rad4::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*11728586107886602578331985848725315229022887949105249281174083471329465317158 + 11728586107886602578331985848725315229022887949105249281174083471329465317158";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_2_921_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp921::CGLDim1Rad2::new();
    let hash = cgl.hash(MSG.to_vec());
    println!("Rust:     {}", hash);

    let expected: &str = "i*39230057073311561979904713106661291660670794873352665774512435740425423658735 + 28281974817376390026165203588145415349405195023788259912005133424119065028093";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_921_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl = thp921::CGLDim1Rad4::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*39695488926984588049869419990040957434879773680406733879881117127139183500293 + 39695488926984588049869419990040957434879773680406733879881117127139183500293";
    println!("SageMath: {}", expected);
}

fn dimension_two_rad_2_127_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp127::CGLDim2Rad2::new();
    let (h1, h2, h3) = cgl.hash(MSG.to_vec(), 3);

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let exp1 = "i*51270636009212854362226904376444348492 + 118854541911646677142182896694474475917";
    let exp2 = "i*37435301830601681338962451584023412917 + 153456715888985573907459419269142792661";
    let exp3 = "i*103658172645780440076386144196752149867 + 70551901995418187699403785643888592714";

    println!("Sage hash: {}", exp1);
    println!("           {0}\n           {1}", exp2, exp3);
}

fn main() {
    println!("================================================================================");
    println!("                  Dimension One CGL with p = 79*2^247 - 1");
    println!("================================================================================");

    dimension_one_rad_2_254_example();
    println!();
    dimension_one_rad_4_254_example();
    println!("\n");

    println!("================================================================================");
    println!("                  Dimension One CGL with p = 2^255 - 921");
    println!("================================================================================");

    dimension_one_rad_2_921_example();
    println!();
    dimension_one_rad_4_921_example();
    println!("\n");

    println!("================================================================================");
    println!("                  Dimension Two CGL with p = 2^127 - 1");
    println!("================================================================================");

    dimension_two_rad_2_127_example();
    println!("\n");
}
