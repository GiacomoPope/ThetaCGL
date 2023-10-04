#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub mod Fp127 {
    const N: usize = 2;
    const BITLEN: usize = 127;
    const MODULUS: [u64; N] = [0xFFFFFFFFFFFFFFFF, 0x6BFFFFFFFFFFFFFF];
    const HALF_MODULUS: [u64; N] = [0x0000000000000000, 0x3600000000000000];
    const R_VAL: [u64; N] = [0x0000000000000002, 0x2800000000000000];
    const MINUS_R_VAL: [u64; N] = [0xFFFFFFFFFFFFFFFD, 0x43FFFFFFFFFFFFFF];
    const DR_VAL: [u64; N] = [0x0000000000000004, 0x5000000000000000];
    const TR_VAL: [u64; N] = [0x0000000000000007, 0x0C00000000000000];
    const QR_VAL: [u64; N] = [0x0000000000000009, 0x3400000000000000];
    const R2_VAL: [u64; N] = [0x425ED097B425ED0F, 0x0AD097B425ED097B];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [0x936CB7E82C388EF5, 0x4A4F426C22137BB3];
    const TDEC_VAL: [u64; N] = [0x5ED097B425ED097B, 0x1C00000000000002];
    const SQRT_EH: [u8; 1] = [27];
    const SQRT_EL: usize = 24;
    const FOURTH_ROOT_EH: [u8; 2] = [16, 13];
    const FOURTH_ROOT_EL: usize = 23;
    const P1: u64 = 3623878655;
    const P1DIV_M: u64 = 3416063723386606283;

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
    const SQRT_EH: [u8; 2] = [15, 2];
    const SQRT_EL: usize = 49;
    const FOURTH_ROOT_EH: [u8; 3] = [16, 7, 1];
    const FOURTH_ROOT_EL: usize = 48;
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

pub mod Fp127Ext {
    use super::Fp127::Fp;
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
