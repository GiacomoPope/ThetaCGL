#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub mod Fp127 {
    pub use crate::finitefield::gf_127_m64::Gf127;
    pub type Fp = Gf127;
}

pub mod Fp127Old {
    const N: usize = 2;
    const BITLEN: usize = 127;
    const MODULUS: [u64; N] = [0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF];
    const HALF_MODULUS: [u64; N] = [0x0000000000000000, 0x4000000000000000];
    const R_VAL: [u64; N] = [0x0000000000000002, 0x0000000000000000];
    const MINUS_R_VAL: [u64; N] = [0xFFFFFFFFFFFFFFFD, 0x7FFFFFFFFFFFFFFF];
    const DR_VAL: [u64; N] = [0x0000000000000004, 0x0000000000000000];
    const TR_VAL: [u64; N] = [0x0000000000000006, 0x0000000000000000];
    const QR_VAL: [u64; N] = [0x0000000000000008, 0x0000000000000000];
    const R2_VAL: [u64; N] = [0x0000000000000004, 0x0000000000000000];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [0x0000000000000100, 0x0000000000000000];
    const TDEC_VAL: [u64; N] = [0x0000000000000000, 0x0000000000000002];
    const WIN_LEN: usize = 4;
    const SQRT_EH: [u8; 1] = [2];
    const SQRT_EL: usize = 31;
    const FOURTH_ROOT_EH: [u8; 1] = [1];
    const FOURTH_ROOT_EL: usize = 31;
    const P1: u64 = 4294967295;
    const P1DIV_M: u64 = 1;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

// 79*2**247 - 1
pub mod Fp254 {
    const N: usize = 4;
    const BITLEN: usize = 254;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x277FFFFFFFFFFFFF,
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x13C0000000000000,
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000006,
        0x0000000000000000,
        0x0000000000000000,
        0x1300000000000000,
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFF9,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x147FFFFFFFFFFFFF,
    ];
    const DR_VAL: [u64; N] = [
        0x000000000000000C,
        0x0000000000000000,
        0x0000000000000000,
        0x2600000000000000,
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000000013,
        0x0000000000000000,
        0x0000000000000000,
        0x1180000000000000,
    ];
    const QR_VAL: [u64; N] = [
        0x0000000000000019,
        0x0000000000000000,
        0x0000000000000000,
        0x2480000000000000,
    ];
    const R2_VAL: [u64; N] = [
        0xCF6474A8819EC913,
        0x1D2A2067B23A5440,
        0x8819EC8E951033D9,
        0x0CA3A5440CF6474A,
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0xC55EA6D998EC13EE,
        0xF1FCF3BEA10D1E49,
        0xF08483F5FFE96735,
        0x1531789831319D25,
    ];
    const TDEC_VAL: [u64; N] = [
        0x1D2A2067B23A5440,
        0x8819EC8E951033D9,
        0x7B23A5440CF6474A,
        0x2000000000000006,
    ];
    // const WIN_LEN: usize = 5;
    // const SQRT_EH: [u8; 2] = [15, 2];
    // const SQRT_EL: usize = 49;
    // const FOURTH_ROOT_EH: [u8; 3] = [16, 7, 1];
    // const FOURTH_ROOT_EL: usize = 48;
    const WIN_LEN: usize = 4;
    const SQRT_EH: [u8; 2] = [14, 9];
    const SQRT_EL: usize = 61;
    const FOURTH_ROOT_EH: [u8; 2] = [15, 4];
    const FOURTH_ROOT_EL: usize = 61;

    const P1: u64 = 2650800127;
    const P1DIV_M: u64 = 11441651398765969958;

    crate::finitefield::fp_gen::define_fp_core! {}

    #[cfg(test)]
    mod tests {
        crate::finitefield::fp_gen::define_fp_tests! {}
    }
}

pub mod Fp921 {
    pub use crate::finitefield::gf255_m64::GF255;
    pub type Fp = GF255<921>;
}

pub mod Fp5248 {
    pub use crate::finitefield::gf5_248_m64::GF5_248;
    pub type Fp = GF5_248;
}

pub mod Fp64 {
    pub use crate::finitefield::gf64_257::GFp;
    pub type Fp = GFp;
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

pub mod Fp127OldExt {
    use super::Fp127Old::Fp;
    const NQR_RE: Fp = Fp::new([0xE1F195E63EF9998E, 0x025681E8E25E14E8]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp254Ext {
    use super::Fp254::Fp;
    const NQR_RE: Fp = Fp::new([
        0xCA918CFCEB423A35,
        0xF72CEA622F5F8006,
        0x6103129C769699BD,
        0x182D25CEB1AFA29B,
    ]);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}

pub mod Fp921Ext {
    use super::Fp921::Fp;
    const NQR_RE: Fp = Fp::w64le(2, 0, 0, 0);

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

pub mod Fp64Ext {
    use super::Fp64::Fp;

    const NQR_RE: Fp = Fp::from_u64_reduce(0xF3139B9D0738D2D1);

    crate::finitefield::fp2_gen::define_fp2_core! {}
    #[cfg(test)]
    mod tests {
        crate::finitefield::fp2_gen::define_fp2_tests! {}
    }
}
