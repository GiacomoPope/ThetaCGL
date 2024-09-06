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
    let cgl = thp5248::CGLDim1Rad2::new();
    let hash = cgl.hash(MSG.to_vec());
    println!("Rust:     {}", hash);

    let expected: &str = "i*1048790104549868760381845710182137722666314665670240273769600751651562439610 + 516135818785057378604189666983524574774965518529932074051373940387715252763";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_4_5248_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl = thp5248::CGLDim1Rad4::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*1591777228698064412432158294069564584505508595511695847665789421532577214324 + 1680247207835055787499001411562698045894861713928280006839599027958720694872";
    println!("SageMath: {}", expected);
}

fn dimension_one_rad_8_5248_example() {
    println!("Computing using 8-radical isogenies...");
    let cgl = thp5248::CGLDim1Rad8::new();
    let hash = cgl.hash(MSG.to_vec());

    println!("Rust:     {}", hash);

    let expected = "i*623953899702050406121574801940143023523525779863632387363046409469719787819 + 141743768810756281317922111945465696489535516588307365005176069147968563672";
    println!("SageMath: {}", expected);
}

fn dimension_two_rad_2_127_example() {
    println!("Computing using 2-radical isogenies...");
    let cgl = thp127::CGLDim2Rad2::new();
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1 = "i*12320047962523423279486361112665882772 + 108134968326043731154961119248302170639";
    let ex2 = "i*15803407118249139560079599117651420281 + 66205839636749269947956199724487136674";
    let ex3 = "i*90655641865476592982577206231671450841 + 127570496827730619947759641012315192423";

    println!("Sage hash: {}", ex1);
    println!("           {0}\n           {1}", ex2, ex3);
}

fn dimension_two_rad_4_127_example() {
    println!("Computing using 4-radical isogenies...");
    let cgl: thp127::CGLDim2Rad4 = thp127::CGLDim2Rad4::new();
    let (h1, h2, h3) = cgl.hash(MSG.to_vec());

    println!("Rust hash: {}", h1);
    println!("           {0}\n           {1}", h2, h3);

    let ex1 = "i*47375232055273767260534253219319124145 + 73346381814437457025019178306100270373";
    let ex2 = "i*100355035873619011046902920138034070025 + 53620892648164048340211136423015102743";
    let ex3 = "i*30161281981729506452828974222531171272 + 38825259926942117650301850823946812854";

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

    let ex1 = "i*9131022502302392460 + 8047160213627555690";
    let ex2 = "i*4009428002838015941 + 8894398257467059120";
    let ex3 = "i*2554299219589503533 + 1205992271966875535";
    let ex4 = "i*13544911288836066678 + 14775662367645206447";
    let ex5 = "i*13250694430739928298 + 10506884765976621283";
    let ex6 = "i*17561772940977990423 + 12895187849928915004";
    let ex7 = "i*17712918206892277537 + 63058281213942031";

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
