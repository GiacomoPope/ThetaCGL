#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

#[allow(unused_macros)]
macro_rules! static_assert {
    ($condition:expr) => {
        let _ = &[()][1 - ($condition) as usize];
    };
}

pub mod fields;
pub mod finitefield;

mod dimension_one;
mod dimension_three;
mod dimension_two;
mod util;

// ============================================================

pub mod thp64 {
    pub type Fp = crate::fields::Fp64::Fp;
    pub type Fq = crate::fields::Fp64Ext::Fp2;

    // Pseudorandom domain computed from sage code
    const a0_re: Fp = Fp::from_u64_reduce(1);
    const a0_im: Fp = Fp::from_u64_reduce(0);

    const a1_re: Fp = Fp::from_u64_reduce(0xDC5596E1A8BD61A6);
    const a1_im: Fp = Fp::from_u64_reduce(0x27975EF03960B799);

    const a2_re: Fp = Fp::from_u64_reduce(0x4C78D44019CAA57B);
    const a2_im: Fp = Fp::from_u64_reduce(0x631CC21ECF776EA1);

    const a3_re: Fp = Fp::from_u64_reduce(0xFC257B234A441A12);
    const a3_im: Fp = Fp::from_u64_reduce(0xCC9BD67E3708856A);

    const a4_re: Fp = Fp::from_u64_reduce(0x2C35E4555BF6F161);
    const a4_im: Fp = Fp::from_u64_reduce(0xA73405376A64ECF8);

    const a5_re: Fp = Fp::from_u64_reduce(0xB1122F78FAA9BFF3);
    const a5_im: Fp = Fp::from_u64_reduce(0xA0F357B1C4183714);

    const a6_re: Fp = Fp::from_u64_reduce(0x009D1194B44A6CE4);
    const a6_im: Fp = Fp::from_u64_reduce(0xBD8B7727A38AE735);

    const a7_re: Fp = Fp::from_u64_reduce(0x3C8E925EE4A4A528);
    const a7_im: Fp = Fp::from_u64_reduce(0x6C671FF6C8B712D2);

    const A0: Fq = Fq::new(&a0_re, &a0_im);
    const A1: Fq = Fq::new(&a1_re, &a1_im);
    const A2: Fq = Fq::new(&a2_re, &a2_im);
    const A3: Fq = Fq::new(&a3_re, &a3_im);
    const A4: Fq = Fq::new(&a4_re, &a4_im);
    const A5: Fq = Fq::new(&a5_re, &a5_im);
    const A6: Fq = Fq::new(&a6_re, &a6_im);
    const A7: Fq = Fq::new(&a7_re, &a7_im);

    crate::dimension_three::define_dim_three_theta_core! {}
}

pub mod thp127 {
    pub type Fp = crate::fields::Fp127::Fp;
    pub type Fq = crate::fields::Fp127Ext::Fp2;

    // Pseudorandom domain computed from sage code
    const X0_re: Fp = Fp::w64le(1, 0);
    const X0_im: Fp = Fp::w64le(0, 0);
    const Z0_re: Fp = Fp::w64le(0xB554F0A6CD56898D, 0x596A8D585A516637);
    const Z0_im: Fp = Fp::w64le(0x55C62A13DD1A684C, 0x26925D8555F401D3);
    const U0_re: Fp = Fp::w64le(0xF9117AFF78A5C9D5, 0x7372B03CEF530598);
    const U0_im: Fp = Fp::w64le(0x160218406DE240B5, 0x1C29C6E16DD093B8);
    const V0_re: Fp = Fp::w64le(0xEBEE280EE44DC74A, 0x3513CD545EE73127);
    const V0_im: Fp = Fp::w64le(0x54191B88482D8D6B, 0x4DFBD67E82B0296E);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);
    const U0: Fq = Fq::new(&U0_re, &U0_im);
    const V0: Fq = Fq::new(&V0_re, &V0_im);

    crate::dimension_two::define_dim_two_theta_core! {}
}

pub mod thp5248 {
    pub type Fp = crate::fields::Fp5248::Fp;
    pub type Fq = crate::fields::Fp5248Ext::Fp2;

    // Theta null point for domain
    const X0_re: Fp = Fp::w64le(1, 0, 0, 0);
    const X0_im: Fp = Fp::w64le(0, 0, 0, 0);
    const Z0_re: Fp = Fp::w64le(
        0xC8989D9AE0732F74,
        0x57A518734EB6E287,
        0x86D2007E41B557BE,
        0x0384E9FADF42A057,
    );
    const Z0_im: Fp = Fp::w64le(
        0x608090B2A570F5A7,
        0x7F97E988D436D484,
        0x67BE08A170FBD685,
        0x014903873860E447,
    );
    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    // Torsion point for 8-radical
    const PX0_re: Fp = Fp::w64le(
        0xD4BA7502003DDF70,
        0x7C56A339012F7167,
        0xFE5CD953AF85BACD,
        0x01AFB6E52F6BF65B,
    );
    const PX0_im: Fp = Fp::w64le(
        0x34319D65B6155F58,
        0xED7E75228452961C,
        0xCC354C142D9CB7DF,
        0x0363BB6DF4CC1A44,
    );
    const PZ0_re: Fp = Fp::w64le(
        0xE61B297101B4CE24,
        0xA1C86E256D899ADA,
        0x353BE1B19CA0F309,
        0x0471CE85213D0B38,
    );
    const PZ0_im: Fp = Fp::w64le(
        0x771D8B45C979B97B,
        0x85DF7E30CC0E1308,
        0x300F98F48A880ED8,
        0x01BCBAF4B16E5B7C,
    );
    const PX0: Fq = Fq::new(&PX0_re, &PX0_im);
    const PZ0: Fq = Fq::new(&PZ0_re, &PZ0_im);

    // sqrt 2 for 8-radical
    const fp_sqrt_2_re: Fp = Fp::w64le(
        0xFF805D2A0D52E912,
        0xED25DC2169473610,
        0xE2973DF03F968969,
        0x013A0F3E1D7C72C5,
    );
    const fp2_sqrt_2: Fq = Fq::new(&fp_sqrt_2_re, &Fp::ZERO);

    // eighth root of unity for 8-radical
    const zeta_8_re: Fp = Fp::w64le(
        0x803FD16AF9568B76,
        0x096D11EF4B5C64F7,
        0x0EB46107E034BB4B,
        0x0462F860F141C69D,
    );
    const zeta_8_im: Fp = Fp::w64le(
        0x803FD16AF9568B76,
        0x096D11EF4B5C64F7,
        0x0EB46107E034BB4B,
        0x0462F860F141C69D,
    );
    const fp2_zeta_8: Fq = Fq::new(&zeta_8_re, &zeta_8_im);

    crate::dimension_one::define_dim_one_theta_core! {}
}

#[cfg(test)]
mod cgl_tests {
    use super::*;
    static MSG: [u8; 256] = [
        1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0,
        0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
        0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0,
        1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1,
        0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1,
    ];

    #[test]
    fn test_dim_one_rad_two() {
        let block_size = 324;
        let cgl = thp5248::CGLDim1Rad2::new(block_size);
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected: &str = "i*1202388955281156761833005543460283826261893173763289830513707326666207689090 + 1017138393049698567064566102452423541172253177172267716238860301237229426197";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_four() {
        let block_size = 324;
        let cgl = thp5248::CGLDim1Rad4::new(block_size);
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected = "i*581317035081276154367619389388237654364317674072512259556056909695832618407 + 1591777228698064412432158294069564584505508595511695847665789421532577214324";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_eight() {
        let block_size = 324;
        let cgl = thp5248::CGLDim1Rad8::new(block_size);
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected = "i*1494595583073012467876291585018799554575832309229632536072774587258798427789 + 961562885477037459496521684366018410270005839775805972818305096331467724997";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_two_rad_two() {
        let block_size = 324;
        let cgl = thp127::CGLDim2Rad2::new(block_size);
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 = "i*20729702215786162559161654607312901632 + 24824195667478872002224501962719953929";
        let ex2 = "i*26788238508836083525728620758542815739 + 79801762727444263675497163220240127218";
        let ex3 = "i*9208197243585817955268805680162455315 + 162984140915487797143139789109763110655";
        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_two_rad_four() {
        let block_size = 324;
        let cgl = thp127::CGLDim2Rad4::new(block_size);
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 = "i*71334330726735141907711867275106993899 + 28253928546883897112004727297228358195";
        let ex2 = "i*99662282611540068483071581849432213565 + 123905390209433693591331379456627201789";
        let ex3 = "i*31562779578250128665683425681981064250 + 21432070052199729518009820487250206136";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_three_rad_two() {
        let block_size = 324;
        let cgl = thp64::CGLDim3Rad2::new(block_size);
        let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 = "i*1326176967233534206 + 11469748429841971151";
        let ex2 = "i*15933490042872426140 + 4344165148600749617";
        let ex3 = "i*7149887378893630726 + 11252021391754383167";
        let ex4 = "i*1991079645619143973 + 13082157033989801270";
        let ex5 = "i*1233115385892712047 + 78263088787484906";
        let ex6 = "i*9848039998401728707 + 14416038328361622008";
        let ex7 = "i*865834649854545355 + 696009373114088875";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
        assert_eq!(ex4, format!("{}", h4));
        assert_eq!(ex5, format!("{}", h5));
        assert_eq!(ex6, format!("{}", h6));
        assert_eq!(ex7, format!("{}", h7));
    }
}
