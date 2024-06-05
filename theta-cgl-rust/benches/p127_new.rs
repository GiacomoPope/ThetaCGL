mod fpcore;
use theta_cgl_rust::fields::Fp127New::Fp;
use theta_cgl_rust::fields::Fp127NewExt::Fp2;

const FP_NAME: &str = "f*2^127 - 1 (new)";
fpcore::define_fp_bench! {}
