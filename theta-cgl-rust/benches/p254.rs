mod fpcore;
use theta_cgl_rust::fields::Fp254::Fp;
use theta_cgl_rust::fields::Fp254Ext::Fp2;

const FP_NAME: &str = "79*2^247 - 1";
fpcore::define_fp_bench!{}