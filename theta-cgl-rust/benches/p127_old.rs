mod fpcore;
use theta_cgl_rust::fields::Fp127::Fp;
use theta_cgl_rust::fields::Fp127Ext::Fp2;

const FP_NAME: &str = "f*2^127 - 1 (old)";
fpcore::define_fp_bench! {}
