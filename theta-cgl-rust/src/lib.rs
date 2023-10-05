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
mod util;

// ============================================================

pub mod thp127 {
    pub type Fp = crate::fields::Fp127::Fp;
    pub type Fq = crate::fields::Fp127Ext::Fp2;

    //
    // WARNING!!
    // The theta coordinates constructed here are assumed to already be
    // in Montgomery form!!!
    //
    const X0_re: Fp = Fp::new([0xFFFFFFFFFFFFFFFB, 0x1BFFFFFFFFFFFFFF]);
    const X0_im: Fp = Fp::new([0, 0]);
    const Z0_re: Fp = Fp::new([0x9B0FBD70A4230C25, 0x00966576AE856B61]);
    const Z0_im: Fp = Fp::new([0x9B0FBD70A4230C25, 0x00966576AE856B61]);
    const U0_re: Fp = Fp::new([0x9B0FBD70A4230C25, 0x00966576AE856B61]);
    const U0_im: Fp = Fp::new([0x9B0FBD70A4230C25, 0x00966576AE856B61]);
    const V0_re: Fp = Fp::new([0, 0]);
    const V0_im: Fp = Fp::new([0xFFFFFFFFFFFFFFFB, 0x1BFFFFFFFFFFFFFF]);

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
    const X0_im: Fp = Fp::new([0x220CDBD2842A931B, 0xD2C4D7CF8280D5F4, 0x6E0784A5D4A77E8D, 0x056433FC148D3399]);
    const Z0_re: Fp = Fp::new([0xFFFFFFFFFFFFFFF9, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x147FFFFFFFFFFFFF]);
    const Z0_im: Fp = Fp::new([0x0000000000000006, 0x0000000000000000, 0x0000000000000000, 0x1300000000000000]);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    crate::dimension_one::define_dim_one_theta_core! {}
}

pub mod thp921 {
    pub type Fp = crate::fields::Fp921::Fp;
    pub type Fq = crate::fields::Fp921Ext::Fp2;

    const X0_re: Fp = Fp::w64le(0, 0, 0, 0);
    const X0_im: Fp = Fp::w64le(0xFF2EAAF0CFC9C0C4, 0x31418C3335B9840B, 0xFC83FEE1060A26D4, 0x507A4C86598511FB);
    const Z0_re: Fp = Fp::w64le(0xFFFFFFFFFFFFFC66, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF);
    const Z0_im: Fp = Fp::w64le(1, 0, 0, 0);

    const X0: Fq = Fq::new(&X0_re, &X0_im);
    const Z0: Fq = Fq::new(&Z0_re, &Z0_im);

    crate::dimension_one::define_dim_one_theta_core! {}
}
