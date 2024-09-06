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
        let cgl = thp5248::CGLDim1Rad2::new();
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected: &str = "i*1048790104549868760381845710182137722666314665670240273769600751651562439610 + 516135818785057378604189666983524574774965518529932074051373940387715252763";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_four() {
        let cgl = thp5248::CGLDim1Rad4::new();
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected = "i*1591777228698064412432158294069564584505508595511695847665789421532577214324 + 1680247207835055787499001411562698045894861713928280006839599027958720694872";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_eight() {
        let cgl = thp5248::CGLDim1Rad8::new();
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected = "i*623953899702050406121574801940143023523525779863632387363046409469719787819 + 141743768810756281317922111945465696489535516588307365005176069147968563672";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_two_rad_two() {
        let cgl = thp127::CGLDim2Rad2::new();
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 =
            "i*12320047962523423279486361112665882772 + 108134968326043731154961119248302170639";
        let ex2 =
            "i*15803407118249139560079599117651420281 + 66205839636749269947956199724487136674";
        let ex3 =
            "i*90655641865476592982577206231671450841 + 127570496827730619947759641012315192423";
        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_two_rad_four() {
        let cgl = thp127::CGLDim2Rad4::new();
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 =
            "i*47375232055273767260534253219319124145 + 73346381814437457025019178306100270373";
        let ex2 =
            "i*100355035873619011046902920138034070025 + 53620892648164048340211136423015102743";
        let ex3 =
            "i*30161281981729506452828974222531171272 + 38825259926942117650301850823946812854";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_three_rad_two() {
        let cgl = thp64::CGLDim3Rad2::new();
        let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 = "i*9131022502302392460 + 8047160213627555690";
        let ex2 = "i*4009428002838015941 + 8894398257467059120";
        let ex3 = "i*2554299219589503533 + 1205992271966875535";
        let ex4 = "i*13544911288836066678 + 14775662367645206447";
        let ex5 = "i*13250694430739928298 + 10506884765976621283";
        let ex6 = "i*17561772940977990423 + 12895187849928915004";
        let ex7 = "i*17712918206892277537 + 63058281213942031";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
        assert_eq!(ex4, format!("{}", h4));
        assert_eq!(ex5, format!("{}", h5));
        assert_eq!(ex6, format!("{}", h6));
        assert_eq!(ex7, format!("{}", h7));
    }
}
