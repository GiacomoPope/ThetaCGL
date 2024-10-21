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

    const a1_re: Fp = Fp::from_u64_reduce(0xD9B9E6A312EB10E3);
    const a1_im: Fp = Fp::from_u64_reduce(0x0FACB7058EB3B138);

    const a2_re: Fp = Fp::from_u64_reduce(0x2D68C7DE2DDE4AB0);
    const a2_im: Fp = Fp::from_u64_reduce(0x7C628F1A55991804);

    const a3_re: Fp = Fp::from_u64_reduce(0x00A84F95E8123329);
    const a3_im: Fp = Fp::from_u64_reduce(0xC39E1A6A8C5F5C81);

    const a4_re: Fp = Fp::from_u64_reduce(0x7A562B50E2981860);
    const a4_im: Fp = Fp::from_u64_reduce(0x17F3AD3F84EED1DA);

    const a5_re: Fp = Fp::from_u64_reduce(0x60830A6D51CE3888);
    const a5_im: Fp = Fp::from_u64_reduce(0xD160A97019010209);

    const a6_re: Fp = Fp::from_u64_reduce(0xFFCF8077B8C8B776);
    const a6_im: Fp = Fp::from_u64_reduce(0x87FF9C2594B7F822);

    const a7_re: Fp = Fp::from_u64_reduce(0xE593E5EA1DD3AF91);
    const a7_im: Fp = Fp::from_u64_reduce(0xE9D550F205F88961);

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
    const Z0_re: Fp = Fp::w64le(0x325FDC555723D5A0, 0x77F5FC2B726DB731);
    const Z0_im: Fp = Fp::w64le(0xBA99785C7BFAEE3A, 0x720B2F8B417528C9);
    const U0_re: Fp = Fp::w64le(0xF8ADFD2EF3323E66, 0x6EE9E8B0A9787DCC);
    const U0_im: Fp = Fp::w64le(0x92B3EC5423248820, 0x1E00DE68A0A2F235);
    const V0_re: Fp = Fp::w64le(0x1BBB14606A71BFFA, 0x2CAADC55E9E1B664);
    const V0_im: Fp = Fp::w64le(0x553F518A7B0E16C7, 0x242A69FCE917627D);

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
        0xE182300064425B3B,
        0xD3D095D71A78C89D,
        0xFBCAA24A766E05BD,
        0x016A0FE894BE77D4,
    );
    const Z0_im: Fp = Fp::w64le(
        0xC4846452776BAC87,
        0x3360EA9B0E196E3C,
        0x7F8FCF9B7D2D4644,
        0x005529F061E8B0F2,
    );
    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    // Torsion point for 8-radical
    const PX0_re: Fp = Fp::w64le(
        0x011265EF861ACDC0,
        0x1927B02D4B32DFB8,
        0x2FD7CB1B8B6BD634,
        0x033F0B4B3D1715CB,
    );
    const PX0_im: Fp = Fp::w64le(
        0xA101F9797BC3EFC7,
        0x2BC8C6892109FE08,
        0x6DE30819F6F9A5B9,
        0x04AA90FBC69C31A6,
    );
    const PZ0_re: Fp = Fp::w64le(
        0x95A19B58F91D0CEC,
        0x2FFEFC7FDE2E8354,
        0xAFDE75D34AF96FBF,
        0x03E4244162713ECB,
    );
    const PZ0_im: Fp = Fp::w64le(
        0xA5B0DEF086D16794,
        0x64ED1CB692F5C089,
        0xD669E53DC5A120E3,
        0x021FD07F173A4C97,
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

        let expected: &str = "i*1175115321110588234636748827639696740191085662044029690205157067136297489851 + 1733123938962694390472365711257600165904333930518679618652142301506945904479";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_four() {
        let block_size = 324;
        let cgl = thp5248::CGLDim1Rad4::new(block_size);
        let hash = cgl.hash(MSG.to_vec());

        let expected: &str = "i*639708434017176960680363788011691951478630737464825183378028357205845236539 + 2111092790111734845126313069227408483857723867262697296810704139886318080888";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_eight() {
        let block_size = 324;
        let cgl = thp5248::CGLDim1Rad8::new(block_size);
        let hash = cgl.hash(MSG.to_vec());

        let expected: &str = "i*1713152820722042121560483199973967855022415860828735853279579284772517053665 + 1614205432150870091530419952171751517745377136920762256637389075657697131180";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_two_rad_two() {
        let block_size = 324;
        let cgl = thp127::CGLDim2Rad2::new(block_size);
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        let ex1: &str =
            "i*55352596873554344554080451593955087743 + 782814348061583778580831952510330245";
        let ex2: &str =
            "i*3636415493409036562625612378128461539 + 93312549852846027507575952345482069462";
        let ex3: &str =
            "i*116984420811334091568360555719211342140 + 115322330464924214960091192691425793045";
        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_two_rad_four() {
        let block_size = 324;
        let cgl = thp127::CGLDim2Rad4::new(block_size);
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        let ex1: &str =
            "i*24343354270362226702211756234962239688 + 165218543846196481307174152989384967825";
        let ex2: &str =
            "i*106957632277458656937006174938906787144 + 148088445431981206886473758973254407305";
        let ex3: &str =
            "i*79492710385759828854917783756632597417 + 22800526872756505643884305741401170041";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_three_rad_two() {
        let block_size = 324;
        let cgl = thp64::CGLDim3Rad2::new(block_size);
        let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

        let ex1: &str = "i*14715636858552704983 + 15452868413087287111";
        let ex2: &str = "i*420104897209182168 + 14547627927149152291";
        let ex3: &str = "i*7553490361471112083 + 16421010721108912633";
        let ex4: &str = "i*8376455609674458376 + 13944539442861068611";
        let ex5: &str = "i*7636771789307779928 + 11852614693545011447";
        let ex6: &str = "i*11092513195598203570 + 6788899222922492439";
        let ex7: &str = "i*10879969565550905398 + 14423116945163197625";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
        assert_eq!(ex4, format!("{}", h4));
        assert_eq!(ex5, format!("{}", h5));
        assert_eq!(ex6, format!("{}", h6));
        assert_eq!(ex7, format!("{}", h7));
    }
}
