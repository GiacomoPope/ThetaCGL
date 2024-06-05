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
mod dimension_two;
mod dimension_three;
mod util;

// ============================================================

pub mod thp64 {
    pub type Fp = crate::fields::Fp64::Fp;
    pub type Fq = crate::fields::Fp64Ext::Fp2;

    //
    // WARNING!!
    // The theta coordinates constructed here are assumed to already be
    // in Montgomery form!!!
    //
    const X0_re: Fp = Fp::from_u64_reduce(0);
    const X0_im: Fp = Fp::from_u64_reduce(0x748E28D491B95FB6);

    const Z0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const Z0_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);

    const U0_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const U0_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);

    const V0_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const V0_im: Fp = Fp::from_u64_reduce(0);

    const G_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const G_im: Fp = Fp::from_u64_reduce(0xFFFFFFFFFFFFFEFD);

    const H_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const H_im: Fp = Fp::from_u64_reduce(0);

    const I_re: Fp = Fp::from_u64_reduce(0x8B71D72B6E469F49);
    const I_im: Fp = Fp::from_u64_reduce(0);

    const J_re: Fp = Fp::from_u64_reduce(0x0000000000000002);
    const J_im: Fp = Fp::from_u64_reduce(0x0000000000000002);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);
    const U0: Fq = Fq::new(&U0_re, &U0_im);
    const V0: Fq = Fq::new(&V0_re, &V0_im);

    const G0: Fq = Fq::new(&G_re, &G_im);
    const H0: Fq = Fq::new(&H_re, &H_im);
    const I0: Fq = Fq::new(&I_re, &I_im);
    const J0: Fq = Fq::new(&J_re, &J_im);

    crate::dimension_three::define_dim_three_theta_core! {}
}

pub mod thp127 {
    pub type Fp = crate::fields::Fp127::Fp;
    pub type Fq = crate::fields::Fp127Ext::Fp2;

    //
    // WARNING!!
    // The theta coordinates constructed here are assumed to already be
    // in Montgomery form!!!
    //
    const X0_re: Fp = Fp::new([0xFFFFFFFFFFFFFFFB, 0x7FFFFFFFFFFFFFFF]);
    const X0_im: Fp = Fp::new([0, 0]);
    const Z0_re: Fp = Fp::new([0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFD]);
    const Z0_im: Fp = Fp::new([0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFD]);
    const U0_re: Fp = Fp::new([0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFD]);
    const U0_im: Fp = Fp::new([0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFD]);
    const V0_re: Fp = Fp::new([0, 0]);
    const V0_im: Fp = Fp::new([0xFFFFFFFFFFFFFFFB, 0x7FFFFFFFFFFFFFFF]);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);
    const U0: Fq = Fq::new(&U0_re, &U0_im);
    const V0: Fq = Fq::new(&V0_re, &V0_im);

    crate::dimension_two::define_dim_two_theta_core! {}
}

pub mod thp127new {
    pub type Fp = crate::fields::Fp127New::Fp;
    pub type Fq = crate::fields::Fp127NewExt::Fp2;

    //
    // WARNING!!
    // The theta coordinates constructed here are assumed to already be
    // in Montgomery form!!!
    //

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

pub mod thp254 {
    pub type Fp = crate::fields::Fp254::Fp;
    pub type Fq = crate::fields::Fp254Ext::Fp2;

    //
    // WARNING!!
    // The theta coordinates constructed here are assumed to already be
    // in Montgomery form!!!
    //
    const X0_re: Fp = Fp::new([0, 0, 0, 0]);
    const X0_im: Fp = Fp::new([
        0x220CDBD2842A931B,
        0xD2C4D7CF8280D5F4,
        0x6E0784A5D4A77E8D,
        0x056433FC148D3399,
    ]);
    const Z0_re: Fp = Fp::new([
        0xFFFFFFFFFFFFFFF9,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x147FFFFFFFFFFFFF,
    ]);
    const Z0_im: Fp = Fp::new([
        0x0000000000000006,
        0x0000000000000000,
        0x0000000000000000,
        0x1300000000000000,
    ]);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    crate::dimension_one::define_dim_one_theta_core! {}
}

pub mod thp921 {
    pub type Fp = crate::fields::Fp921::Fp;
    pub type Fq = crate::fields::Fp921Ext::Fp2;

    const X0_re: Fp = Fp::w64le(0, 0, 0, 0);
    const X0_im: Fp = Fp::w64le(
        0xFF2EAAF0CFC9C0C4,
        0x31418C3335B9840B,
        0xFC83FEE1060A26D4,
        0x507A4C86598511FB,
    );
    const Z0_re: Fp = Fp::w64le(
        0xFFFFFFFFFFFFFC66,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x7FFFFFFFFFFFFFFF,
    );
    const Z0_im: Fp = Fp::w64le(1, 0, 0, 0);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    crate::dimension_one::define_dim_one_theta_core! {}
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
