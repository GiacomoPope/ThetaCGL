#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp127old;
use theta_cgl_rust::thp254;
use theta_cgl_rust::thp5248;
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

fn dimension_one_rad_2_5248_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp5248::CGLDim1Rad2::new();
    let hash = cgl.hash(MSG.to_vec());
    println!("Rust:     {}", hash);

    let expected: &str = "i*581317035081276154367619389388237654364317674072512259556056909695832618407 + 1591777228698064412432158294069564584505508595511695847665789421532577214324";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_5248_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl = thp5248::CGLDim1Rad4::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*1984116653337552224850207110995407685821969784456245019753456222927051983734 + 1984116653337552224850207110995407685821969784456245019753456222927051983734";
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

fn dimension_two_rad_2_127_old_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp127old::CGLDim2Rad2::new();
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1 = "i*51270636009212854362226904376444348492 + 118854541911646677142182896694474475917";
    let ex2 = "i*37435301830601681338962451584023412917 + 153456715888985573907459419269142792661";
    let ex3 = "i*103658172645780440076386144196752149867 + 70551901995418187699403785643888592714";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_two_rad_4_127_old_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl = thp127old::CGLDim2Rad4::new();
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

    let ex1 = "i*8195311847842562697 + 15367522950815623545";
    let ex2 = "i*5509503292246740166 + 3363830550895998170";
    let ex3 = "i*6222483436527344390 + 5906729632957704994";
    let ex4 = "i*9457393622803510435 + 14560144086626072899";
    let ex5 = "i*12039521041378489879 + 15342103277418457660";
    let ex6 = "i*7403012461811820402 + 12518067020380020935";
    let ex7 = "i*14662549957566082920 + 1196490411222064188";

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
    println!("                  Dimension One CGL with p = 5*2^258 - 1");
    println!("================================================================================");

    dimension_one_rad_2_5248_example();
    println!();
    dimension_one_rad_4_5248_example();
    println!("\n");

    println!("================================================================================");
    println!("                    Dimension Two CGL with p = 2^127 - 1");
    println!("================================================================================");

    dimension_two_rad_2_127_example();
    println!();
    dimension_two_rad_4_127_example();
    println!("\n");

    println!("================================================================================");
    println!("                 Dimension Two CGL with p = 2^127 - 1 (old)");
    println!("================================================================================");

    dimension_two_rad_2_127_old_example();
    println!();
    dimension_two_rad_4_127_old_example();
    println!("\n");

    println!("================================================================================");
    println!("                  Dimension Three CGL with p = 2^64 - 257");
    println!("================================================================================");

    dimension_three_rad_2_64_example();
    println!("\n");
}
