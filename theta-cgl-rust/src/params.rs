#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub mod p127 {
    const BITLEN: usize = 127;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000, 0x4000000000000000
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000002, 0x0000000000000000
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFFD, 0x7FFFFFFFFFFFFFFF
    ];
    const DR_VAL: [u64; N] = [
        0x0000000000000004, 0x0000000000000000
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000000006, 0x0000000000000000
    ];
    const QR_VAL: [u64; N] = [
        0x0000000000000008, 0x0000000000000000
    ];
    const R2_VAL: [u64; N] = [
        0x0000000000000004, 0x0000000000000000
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0x0000000000000100, 0x0000000000000000
    ];
    const TDEC_VAL: [u64; N] = [
        0x0000000000000000, 0x0000000000000002
    ];
    const SQRT_EH: [u8; 1] = [
        1
    ];
    const SQRT_EL: usize = 25;
    const P1: u64 = 4294967295;
    const P1DIV_M: u64 = 4294967297;
    const NQR_RE_VAL: [u64; N] = [
        0xB5A7C8E687BA167A, 0x349CD5EFC40B4A42
    ];

    crate::fpcore::define_fp_core!{}

    #[cfg(test)]
    mod tests {
        crate::fpcore::define_fp_tests!{}
    }
}

// 79*2**247 - 1
pub mod p254 {
    const BITLEN: usize = 254;
    const MODULUS: [u64; N] = [
        0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x277FFFFFFFFFFFFF
    ];
    const HALF_MODULUS: [u64; N] = [
        0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x13C0000000000000
    ];
    const R_VAL: [u64; N] = [
        0x0000000000000006, 0x0000000000000000, 0x0000000000000000, 0x1300000000000000
    ];
    const MINUS_R_VAL: [u64; N] = [
        0xFFFFFFFFFFFFFFF9, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x147FFFFFFFFFFFFF
    ];
    const DR_VAL: [u64; N] = [
        0x000000000000000C, 0x0000000000000000, 0x0000000000000000, 0x2600000000000000
    ];
    const TR_VAL: [u64; N] = [
        0x0000000000000013, 0x0000000000000000, 0x0000000000000000, 0x1180000000000000
    ];
    const QR_VAL: [u64; N] = [
        0x0000000000000019, 0x0000000000000000, 0x0000000000000000, 0x2480000000000000
    ];
    const R2_VAL: [u64; N] = [
        0xCF6474A8819EC913, 0x1D2A2067B23A5440, 0x8819EC8E951033D9, 0x0CA3A5440CF6474A
    ];
    const P0I: u64 = 1;
    const TFIXDIV_VAL: [u64; N] = [
        0xC55EA6D998EC13EE, 0xF1FCF3BEA10D1E49, 0xF08483F5FFE96735, 0x1531789831319D25
    ];
    const TDEC_VAL: [u64; N] = [
        0x1D2A2067B23A5440, 0x8819EC8E951033D9, 0x7B23A5440CF6474A, 0x2000000000000006
    ];
    const SQRT_EH: [u8; 2] = [
        15, 2
    ];
    const SQRT_EL: usize = 49;
    const P1: u64 = 2650800127;
    const P1DIV_M: u64 = 11441651398765969958;
    const NQR_RE_VAL: [u64; N] = [
        0xCA918CFCEB423A35, 0xF72CEA622F5F8006, 0x6103129C769699BD, 0x182D25CEB1AFA29B
    ];

    crate::fpcore::define_fp_core!{}

    #[cfg(test)]
    mod tests {
        crate::fpcore::define_fp_tests!{}
    }
}
