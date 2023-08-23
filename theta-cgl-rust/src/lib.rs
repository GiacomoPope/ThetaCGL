#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

mod utils64;
mod fpcore;
mod dimension_one;
mod dimension_two;
mod util;

pub mod params;

// ============================================================

pub mod thp127 {
    pub type Fq = crate::params::p127::Fp2;
    crate::dimension_two::define_dim_two_theta_core!{}
}

pub mod thp254 {
    pub type Fp = crate::params::p254::Fp;
    pub type Fq = crate::params::p254::Fp2;
    crate::dimension_one::define_dim_one_theta_core!{}
}
