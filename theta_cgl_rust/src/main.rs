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

    let expected: &str = "i*1175115321110588234636748827639696740191085662044029690205157067136297489851 + 1733123938962694390472365711257600165904333930518679618652142301506945904479";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_5248_example() {
    println!("Computing using 4-radical isogenies...");
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad4::new(block_size);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected: &str = "i*639708434017176960680363788011691951478630737464825183378028357205845236539 + 2111092790111734845126313069227408483857723867262697296810704139886318080888";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_8_5248_example() {
    println!("Computing using 8-radical isogenies...");
    let block_size = 324;
    let cgl = thp5248::CGLDim1Rad8::new(block_size);
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected: &str = "i*1713152820722042121560483199973967855022415860828735853279579284772517053665 + 1614205432150870091530419952171751517745377136920762256637389075657697131180";
    println!("SageMath: {}", expected);
}

fn dimension_two_rad_2_127_example() {
    println!("Computing using 2-radical isogenies...");
    let block_size = 324;
    let cgl = thp127::CGLDim2Rad2::new(block_size);
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1: &str =
        "i*55352596873554344554080451593955087743 + 782814348061583778580831952510330245";
    let ex2: &str =
        "i*3636415493409036562625612378128461539 + 93312549852846027507575952345482069462";
    let ex3: &str =
        "i*116984420811334091568360555719211342140 + 115322330464924214960091192691425793045";

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

    let ex1: &str =
        "i*24343354270362226702211756234962239688 + 165218543846196481307174152989384967825";
    let ex2: &str =
        "i*106957632277458656937006174938906787144 + 148088445431981206886473758973254407305";
    let ex3: &str =
        "i*79492710385759828854917783756632597417 + 22800526872756505643884305741401170041";

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

    let ex1: &str = "i*14715636858552704983 + 15452868413087287111";
    let ex2: &str = "i*420104897209182168 + 14547627927149152291";
    let ex3: &str = "i*7553490361471112083 + 16421010721108912633";
    let ex4: &str = "i*8376455609674458376 + 13944539442861068611";
    let ex5: &str = "i*7636771789307779928 + 11852614693545011447";
    let ex6: &str = "i*11092513195598203570 + 6788899222922492439";
    let ex7: &str = "i*10879969565550905398 + 14423116945163197625";

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
