mod fpcore;
use theta_cgl_rust::fields::Fp5248::Fp;
use theta_cgl_rust::fields::Fp5248Ext::Fp2;

const FP_NAME: &str = "5*2^248 - 1";
fpcore::define_fp_bench! {}
