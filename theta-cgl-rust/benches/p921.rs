mod fpcore;
use theta_cgl_rust::fields::Fp921::Fp;
use theta_cgl_rust::fields::Fp921Ext::Fp2;

const FP_NAME: &str = "2**255 - 921";
fpcore::define_fp_bench!{}