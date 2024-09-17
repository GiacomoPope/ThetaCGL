#![allow(non_snake_case)]

use theta_cgl_rust::thp127;
use theta_cgl_rust::thp5248;
use theta_cgl_rust::thp64;

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

fn dimension_one_rad_2_5248_example() {
    println!("Computing using 2-radical isogenies...");
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad2::new(block_size);
    let hash = cgl.hash(MSG.to_vec());
    println!("Rust:     {}", hash); 

    let expected: &str = "i*1202388955281156761833005543460283826261893173763289830513707326666207689090 + 1017138393049698567064566102452423541172253177172267716238860301237229426197";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_5248_example() {
    println!("Computing using 4-radical isogenies...");
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad4::new(block_size);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*581317035081276154367619389388237654364317674072512259556056909695832618407 + 1591777228698064412432158294069564584505508595511695847665789421532577214324";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_8_5248_example() {
    println!("Computing using 8-radical isogenies...");
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad8::new(block_size);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*1494595583073012467876291585018799554575832309229632536072774587258798427789 + 961562885477037459496521684366018410270005839775805972818305096331467724997";
    println!("SageMath: {}", expected);
}

fn dimension_two_rad_2_127_example() {
    println!("Computing using 2-radical isogenies...");
    let block_size = 324;
    let cgl = thp127::CGLDim2Rad2::new(block_size);
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3); 

    let ex1 = "i*20729702215786162559161654607312901632 + 24824195667478872002224501962719953929";
    let ex2 = "i*26788238508836083525728620758542815739 + 79801762727444263675497163220240127218";
    let ex3 = "i*9208197243585817955268805680162455315 + 162984140915487797143139789109763110655";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_two_rad_4_127_example() {
    println!("Computing using 4-radical isogenies...");
    let block_size = 324;
    let cgl: thp127::CGLDim2Rad4 = thp127::CGLDim2Rad4::new(block_size);
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3); 

    let ex1 = "i*71334330726735141907711867275106993899 + 28253928546883897112004727297228358195";
    let ex2 = "i*99662282611540068483071581849432213565 + 123905390209433693591331379456627201789";
    let ex3 = "i*31562779578250128665683425681981064250 + 21432070052199729518009820487250206136";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_three_rad_2_64_example() {
    println!("Computing using 2-radical isogenies...");
    let block_size = 324;
    let cgl = thp64::CGLDim3Rad2::new(block_size);
    let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);
    println!("           {0}\n           {1}", h4, h5);
    println!("           {0}\n           {1}", h6, h7);
    println!("");

    let ex1 = "i*1326176967233534206 + 11469748429841971151";
    let ex2 = "i*15933490042872426140 + 4344165148600749617";
    let ex3 = "i*7149887378893630726 + 11252021391754383167";
    let ex4 = "i*1991079645619143973 + 13082157033989801270";
    let ex5 = "i*1233115385892712047 + 78263088787484906";
    let ex6 = "i*9848039998401728707 + 14416038328361622008";
    let ex7 = "i*865834649854545355 + 696009373114088875";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
    println!("           {0}\n           {1}", ex4, ex5);
    println!("           {0}\n           {1}", ex6, ex7);
}

fn main() {
    println!("================================================================================");
    println!("                  Dimension One CGL with p = 5*2^258 - 1");
    println!("================================================================================");

    dimension_one_rad_2_5248_example();
    println!();
    dimension_one_rad_4_5248_example();
    println!();
    dimension_one_rad_8_5248_example();
    println!("\n");

    println!("================================================================================");
    println!("                    Dimension Two CGL with p = 2^127 - 1");
    println!("================================================================================");

    dimension_two_rad_2_127_example();
    println!();
    dimension_two_rad_4_127_example();
    println!("\n");

    println!("================================================================================");
    println!("                  Dimension Three CGL with p = 2^64 - 257");
    println!("================================================================================");

    dimension_three_rad_2_64_example();
    println!("\n");
}
