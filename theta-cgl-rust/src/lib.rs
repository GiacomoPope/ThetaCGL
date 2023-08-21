#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

mod utils64;
mod fpcore;
mod dimension_one;

pub mod params;

// ============================================================

pub mod thp254 {
    pub type Fp = crate::params::p254::Fp2;
    crate::dimension_one::define_dim_one_theta_core!{}
}
