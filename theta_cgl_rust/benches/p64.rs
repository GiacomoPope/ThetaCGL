mod fpcore;
use theta_cgl_rust::fields::Fp64::Fp;
use theta_cgl_rust::fields::Fp64Ext::Fp2;

const FP_NAME: &str = "2^64 - 2^8 - 1";
fpcore::define_fp_bench! {}
