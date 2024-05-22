#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp254;
use theta_cgl_rust::thp64;
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
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1 = "i*51270636009212854362226904376444348492 + 118854541911646677142182896694474475917";
    let ex2 = "i*37435301830601681338962451584023412917 + 153456715888985573907459419269142792661";
    let ex3 = "i*103658172645780440076386144196752149867 + 70551901995418187699403785643888592714";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_two_rad_4_127_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl: thp127::CGLDim2Rad4 = thp127::CGLDim2Rad4::new();
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1 = "i*14255708726834209014572159019899113927 + 106613803422737905865003815327946998256";
    let ex2 = "i*161915183057086815891740274300276651399 + 13519842835831752693363006384794856";
    let ex3 = "i*105249037341253961638806546690278378374 + 165821113257513419160824899234824339511";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_three_rad_2_64_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp64::CGLDim3Rad2::new();
    let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);
    println!("           {0}\n           {1}", h4, h5);
    println!("           {0}\n           {1}", h6, h7);
    println!("");

    let ex1 = "8108166279742740647*i + 18197105568556342780";
    let ex2 = "7220468887645762090*i + 9379785522510560923";
    let ex3 = "14467525117176772143*i + 8786789766168361980";
    let ex4 = "11299964267018974902*i + 10427085028160938082";
    let ex5 = "5064901717020216585*i + 15260026189057691864";
    let ex6 = "1856945183658058105*i + 3943495375381825177";
    let ex7 = "4167688461421429561*i + 2033720393454306340";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
    println!("           {0}\n           {1}", ex4, ex5);
    println!("           {0}\n           {1}", ex6, ex7);
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
    println!();
    dimension_two_rad_4_127_example();

    println!("================================================================================");
    println!("                  Dimension Three CGL with p = 2^64 - 257");
    println!("================================================================================");

    dimension_three_rad_2_64_example();

    println!("\n");
}
