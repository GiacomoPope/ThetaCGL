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

    const X0_re: Fp = Fp::from_u64_reduce(0);
    const X0_im: Fp = Fp::from_u64_reduce(0x748E28D491B95FB6);
    const Z0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const Z0_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);
    const U0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const U0_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);
    const V0_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const V0_im: Fp = Fp::from_u64_reduce(0);
    const G0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const G0_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);
    const H0_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const H0_im: Fp = Fp::from_u64_reduce(0);
    const I0_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const I0_im: Fp = Fp::from_u64_reduce(0);
    const J0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const J0_im: Fp = Fp::from_u64_reduce(0x0000000000000002);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);
    const U0: Fq = Fq::new(&U0_re, &U0_im);
    const V0: Fq = Fq::new(&V0_re, &V0_im);
    const G0: Fq = Fq::new(&G0_re, &G0_im);
    const H0: Fq = Fq::new(&H0_re, &H0_im);
    const I0: Fq = Fq::new(&I0_re, &I0_im);
    const J0: Fq = Fq::new(&J0_re, &J0_im);

    crate::dimension_three::define_dim_three_theta_core! {}
}

pub mod thp127 {
    pub type Fp = crate::fields::Fp127::Fp;
    pub type Fq = crate::fields::Fp127Ext::Fp2;

    const X0_re: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFD, 0x7FFFFFFFFFFFFFFF);
    const X0_im: Fp = Fp::w64le(0, 0);
    const Z0_re: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFE);
    const Z0_im: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFE);
    const U0_re: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFE);
    const U0_im: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFE);
    const V0_re: Fp = Fp::w64le(0, 0);
    const V0_im: Fp = Fp::w64le(0xFFFFFFFFFFFFFFFD, 0x7FFFFFFFFFFFFFFF);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);
    const U0: Fq = Fq::new(&U0_re, &U0_im);
    const V0: Fq = Fq::new(&V0_re, &V0_im);

    crate::dimension_two::define_dim_two_theta_core! {}
}

pub mod thp5248 {
    pub type Fp = crate::fields::Fp5248::Fp;
    pub type Fq = crate::fields::Fp5248Ext::Fp2;

    const X0_re: Fp = Fp::w64le(0, 0, 0, 0);
    const X0_im: Fp = Fp::w64le(
        0xFF805D2A0D52E912,
        0xED25DC2169473610,
        0xE2973DF03F968969,
        0x013A0F3E1D7C72C5,
    );
    const Z0_re: Fp = Fp::w64le(
        0xFFFFFFFFFFFFFFFE,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x04FFFFFFFFFFFFFF,
    );
    const Z0_im: Fp = Fp::w64le(1, 0, 0, 0);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

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
        let expected: &str = "i*581317035081276154367619389388237654364317674072512259556056909695832618407 + 1591777228698064412432158294069564584505508595511695847665789421532577214324";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_one_rad_four() {
        let cgl = thp5248::CGLDim1Rad4::new();
        let hash = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let expected = "i*1984116653337552224850207110995407685821969784456245019753456222927051983734 + 1984116653337552224850207110995407685821969784456245019753456222927051983734";
        assert_eq!(expected, format!("{}", hash));
    }

    #[test]
    fn test_dim_two_rad_two() {
        let cgl = thp127::CGLDim2Rad2::new();
        let (h1, h2, h3) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 =
            "i*51270636009212854362226904376444348492 + 118854541911646677142182896694474475917";
        let ex2 =
            "i*37435301830601681338962451584023412917 + 153456715888985573907459419269142792661";
        let ex3 =
            "i*103658172645780440076386144196752149867 + 70551901995418187699403785643888592714";
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
            "i*14255708726834209014572159019899113927 + 106613803422737905865003815327946998256";
        let ex2 = "i*161915183057086815891740274300276651399 + 13519842835831752693363006384794856";
        let ex3 =
            "i*105249037341253961638806546690278378374 + 165821113257513419160824899234824339511";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
    }

    #[test]
    fn test_dim_three_rad_two() {
        let cgl = thp64::CGLDim3Rad2::new();
        let (h1, h2, h3, h4, h5, h6, h7) = cgl.hash(MSG.to_vec());

        // TODO: don't compare strings! haha
        let ex1 = "i*8195311847842562697 + 15367522950815623545";
        let ex2 = "i*5509503292246740166 + 3363830550895998170";
        let ex3 = "i*6222483436527344390 + 5906729632957704994";
        let ex4 = "i*9457393622803510435 + 14560144086626072899";
        let ex5 = "i*12039521041378489879 + 15342103277418457660";
        let ex6 = "i*7403012461811820402 + 12518067020380020935";
        let ex7 = "i*14662549957566082920 + 1196490411222064188";

        assert_eq!(ex1, format!("{}", h1));
        assert_eq!(ex2, format!("{}", h2));
        assert_eq!(ex3, format!("{}", h3));
        assert_eq!(ex4, format!("{}", h4));
        assert_eq!(ex5, format!("{}", h5));
        assert_eq!(ex6, format!("{}", h6));
        assert_eq!(ex7, format!("{}", h7));
    }
}
