#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

// Base fields GF(p)

pub mod Fp64 {
    pub use crate::finitefield::gf64_257::GFp;
    pub type Fp = GFp;
}

pub mod Fp127 {
    pub use crate::finitefield::gf_127_m64::Gf127;
    pub type Fp = Gf127;
}

pub mod Fp5248 {
    pub use crate::finitefield::gf5_248_m64::GF5_248;
    pub type Fp = GF5_248;
}

// Extension fields GF(p^2) with modulus x^2 + 1

pub mod Fp64Ext {
    use super::Fp64::Fp;

    const NQR_RE: Fp = Fp::from_u64_reduce(5);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp127Ext {
    use super::Fp127::Fp;
    const NQR_RE: Fp = Fp::w64le(2, 0);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp5248Ext {
    use super::Fp5248::Fp;
    const NQR_RE: Fp = Fp::w64le(5, 0, 0, 0);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
