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
    pub type Fq = crate::fields::Fp127Ext::Fp2;
    crate::dimension_two::define_dim_two_theta_core!{}
}

pub mod thp254 {
    pub type Fp = crate::fields::Fp254::Fp;
    pub type Fq = crate::fields::Fp254Ext::Fp2;
    crate::dimension_one::define_dim_one_theta_core!{}
}

pub mod thp921 {
    pub type Fp = crate::fields::Fp921::Fp;
    pub type Fq = crate::fields::Fp921Ext::Fp2;
    crate::dimension_one::define_dim_one_theta_core!{}
}
