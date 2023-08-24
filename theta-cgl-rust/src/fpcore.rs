#![allow(non_snake_case)]

// This code was written by Thomas Pornin, which is part of
// a larger work in progress implementing another project.
//
// Macro expectations:
//   modulus p is prime and larger than 2^64
//   p = 3 mod 4 (required in square roots)
//   defined constants:
//      BITLEN (usize)            modulus length in bits
//      MODULUS ([u64; N])        modulus p (little-endian)
//      HALF_MODULUS ([u64; N])   (p + 1)/2 (little-endian)
//      R_VAL ([u64; N])          2^(64*N) mod p
//      DR_VAL ([u64; N])         2^(64*N+1) mod p
//      R2_VAL ([u64; N])         2^(2*64*N) mod p
//      P0I (u64)                 -1/p mod 2^64
//      TFIXDIV_VAL               corrective factor for division
//      TDEC_VAL                  2^(64*(2*N-1)) mod p
//      SQRT_EH, SQRT_EL          encoding of (p + 1)/4
//      P1 (u64)                  floor(p / 2^(BITLEN - 32))
//      P1DIV_M (u64)             1 + floor((2^32 - P1)*2^64 / P1)
//      NQR_RE_VAL                NQR_RE + i is a non-square in GF(p^2)
//
//   The macro defines the constant N = ceil(BITLEN / 64).
//
// For divisions:
//    let n1 = floor((2*BITLEN - 34) / 31)
//    let n2 = 2*BITLEN - 31*n1 - 2
//    TFIXDIV_VAL = 2^(33*n1 + 64 - n2 + 2*64*N) mod p
//
// For square roots:
//    let e = (p + 1)/4
//    SQRT_EL (usize) is such that e = 0 mod 2^(5*SQRT_EL)
//    SQRT_EH ([u8; ...]) encodes e/2^(5*SQRT_EL) in base 2^5 (little-endian)
//    It is best for performance is SQRT_EL is as big as possible
//
macro_rules! define_fp_core { () => {

    use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
    use core::convert::TryFrom;
    use rand_core::{CryptoRng, RngCore};
    use crate::utils64::{addcarry_u64, subborrow_u64, umull, umull_add, umull_add2, umull_x2, umull_x2_add, sgnw, lzcnt};
    use num_bigint::{BigInt, Sign};
    use std::fmt;

    // Number of limbs
    const N: usize = (BITLEN + 63) >> 6;

    /// A finite field element. Contents are opaque.
    /// All functions are constant-time.
    ///
    /// A field element x is encoded into bytes by using the unsigned
    /// little-endian convention over the unique representant of x in the
    /// [0..(p-1)] range. There is no sign bit.
    #[derive(Clone, Copy, Debug)]
    pub struct Fp([u64; N]);

    impl Fp {

        // IMPLEMENTATION NOTES
        // --------------------
        //
        // Modulus is p. Each element is represented over N limbs, in base
        // 2^64.
        //
        // Let R = 2^(64*N) mod p. A field element x is represented by
        // the integer x*R mod p, in the [0..(p-1)] range. The limbs are in
        // little-endian order (limb 0 is the least significant).
        //
        // Multiplications use Montgomery multiplication: given x and y,
        // the value (x*y)/R mod p is computed. Since our values are in
        // Montgomery representation, what is computed is really
        // (x*R)*(y*R)/R = x*y*R mod p, which is the correct product. The
        // decoding and encoding functions apply the required conversions; the
        // use of Montgomery representation is not visible to other code using
        // this type.

        pub const ZERO: Self = Self([0u64; N]);
        pub const ONE: Self = Self(R_VAL);
        pub const TWO: Self = Self(DR_VAL);
        pub const THREE: Self = Self(TR_VAL);
        pub const FOUR: Self = Self(QR_VAL);
        pub const MINUS_ONE: Self = Self(MINUS_R_VAL);

        const NQR_RE: Self = Self(NQR_RE_VAL);

        /// Modulus bit length.
        pub const BIT_LENGTH: usize = BITLEN;

        /// Encoding length of a field element (in bytes). All elements
        /// always encode into exactly that many bytes. Encoding is
        /// canonical: a given field element has a unique valid encoding,
        /// and the decoding process verifies that this specific encoding
        /// was used.
        pub const ENCODED_LENGTH: usize = (BITLEN + 7) >> 3;

        // R2 = R^2 mod p = Montgomery representation of R
        const R2: Self = Self(R2_VAL);

        // Multiplier for decode_reduce().
        const TDEC: Self = Self(TDEC_VAL);

        // Corrective factor for division.
        const TFIXDIV: Self = Self(TFIXDIV_VAL);

        /// Create an element by converting the provided integer.
        pub fn from_i32(x: i32) -> Self {
            let sx = (x >> 31) as u32;
            let ax = ((x as u32) ^ sx).wrapping_sub(sx);
            let mut r = Self::from_u64(ax as u64);
            r.set_condneg(sx);
            r
        }

        /// Create an element by converting the provided integer.
        pub fn from_i64(x: i64) -> Self {
            let sx = (x >> 63) as u64;
            let ax = ((x as u64) ^ sx).wrapping_sub(sx);
            let mut r = Self::from_u64(ax);
            r.set_condneg(sx as u32);
            r
        }

        /// Create an element by converting the provided integer.
        pub fn from_u32(x: u32) -> Self {
            Self::from_u64(x as u64)
        }

        /// Create an element by converting the provided integer.
        pub fn from_u64(x: u64) -> Self {
            let mut r = Self::ZERO;
            r.0[0] = x;
            r.set_mul(&Self::R2);
            r
        }

        /// Return 0xFFFFFFFF if this value is zero, or 0x00000000 otherwise.
        pub fn iszero(self) -> u32 {
            let mut x = 0;
            for i in 1..N {
                x |= self.0[i];
            }
            (!sgnw(x | x.wrapping_neg())) as u32
        }

        /// Return 0xFFFFFFFF if this value is equal to rhs, or 0x00000000
        /// otherwise.
        pub fn equals(self, rhs: &Self) -> u32 {
            let mut r = 0u64;
            for i in 0..N {
                r |= self.0[i] ^ rhs.0[i];
            }
            (((r | r.wrapping_neg()) >> 63) as u32).wrapping_sub(1)
        }

        /// Add `rhs` to this value.
        fn set_add(&mut self, rhs: &Self) {
            // raw addition.
            let mut cc1 = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], rhs.0[i], cc1);
                self.0[i] = d;
                cc1 = ee;
            }

            // subtract modulus.
            let mut cc2 = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(self.0[i], MODULUS[i], cc2);
                self.0[i] = d;
                cc2 = ee;
            }

            // add back modulus if the result was negative, i.e. cc1 - cc2 < 0.
            let mm = (cc1 as u64).wrapping_sub(cc2 as u64);
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        /// Subtract `rhs` from this value.
        fn set_sub(&mut self, rhs: &Self) {
            // raw subtraction
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(self.0[i], rhs.0[i], cc);
                self.0[i] = d;
                cc = ee;
            }

            // add back modulus if the result was negative
            let mm = (cc as u64).wrapping_neg();
            cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        /// Negate this value.
        pub fn set_neg(&mut self) {
            // subtract from zero
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(0, self.0[i], cc);
                self.0[i] = d;
                cc = ee;
            }

            // add back the modulus if needed (i.e. if input was non-zero)
            let mm = (cc as u64).wrapping_neg();
            cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        // Perform Montgomery reduction (division by R) on this value.
        // Internal note: if self has proper contents (value less than p), then
        // this necessarily yields a properly reduced value. If self is not
        // properly reduced, then the output is in [0..p] inclusive.
        fn set_montyred(&mut self) {
            for _ in 0..N {
                let f = self.0[0].wrapping_mul(P0I);
                let (_, mut cc) = umull_add(f, MODULUS[0], self.0[0]);
                for i in 1..N {
                    let (d, hi) = umull_add2(f, MODULUS[i], self.0[i], cc);
                    self.0[i - 1] = d;
                    cc = hi;
                }
                self.0[N - 1] = cc;
            }
        }

        /// Multiply this value by `rhs`.
        fn set_mul(&mut self, rhs: &Self) {
            let mut t = Self::ZERO;

            // combined muls + reduction
            let mut cch = 0;
            for i in 0..N {
                let f = rhs.0[i];
                let (lo, mut cc1) = umull_add(f, self.0[0], t.0[0]);
                let g = lo.wrapping_mul(P0I);
                let (_, mut cc2) = umull_add(g, MODULUS[0], lo);
                for j in 1..N {
                    let (d, hi1) = umull_add2(f, self.0[j], t.0[j], cc1);
                    cc1 = hi1;
                    let (d, hi2) = umull_add2(g, MODULUS[j], d, cc2);
                    cc2 = hi2;
                    t.0[j - 1] = d;
                }
                let (d, ee) = addcarry_u64(cc1, cc2, cch);
                t.0[N - 1] = d;
                cch = ee;
            }

            // final reduction: subtract modulus if necessary
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(t.0[i], MODULUS[i], cc);
                t.0[i] = d;
                cc = ee;
            }
            let mm = (cch as u64).wrapping_sub(cc as u64);
            cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(t.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        /// Replace this value with its square.
        pub fn set_square(&mut self) {
            // FIXME: this turns out to be slower than set_mul() on x86_64
            // when N >= 23. This is probably due to the more complicated
            // loop bounds. Full unrolling helps, but can only be done at
            // the crate level:
            //   RUSTFLAGS="-C llvm-args=-unroll-threshold=1200"
            // This impacts all the code in the crate, and is thus
            // probably not a very good idea.

            // Compute the square over integers.
            let mut t = [0u64; N << 1];

            // sum_{i<j} a_i*a_j*2^(64*(i+j)) < 2^(64*(2*N-1))
            // -> t[2*N-1] remains at zero
            let f = self.0[0];
            let (d, mut cc) = umull(f, self.0[1]);
            t[1] = d;
            for j in 2..N {
                let (d, hi) = umull_add(f, self.0[j], cc);
                t[j] = d;
                cc = hi;
            }
            t[N] = cc;
            for i in 1..(N - 1) {
                let f = self.0[i];
                let (d, mut cc) = umull_add(f, self.0[i + 1], t[(i << 1) + 1]);
                t[(i << 1) + 1] = d;
                for j in (i + 2)..N {
                    let (d, hi) = umull_add2(f, self.0[j], t[i + j], cc);
                    t[i + j] = d;
                    cc = hi;
                }
                t[i + N] = cc;
            }

            // Double the partial sum.
            // -> t contains sum_{i!=j} a_i*a_j*2^(64*(i+j))
            let mut cc = 0;
            for i in 1..((N << 1) - 1) {
                let w = t[i];
                let ee = w >> 63;
                t[i] = (w << 1) | cc;
                cc = ee;
            }
            t[(N << 1) - 1] = cc;

            // Add the squares a_i*a_i*w^(64*2*i).
            let mut cc = 0;
            for i in 0..N {
                let (lo, hi) = umull(self.0[i], self.0[i]);
                let (d0, ee) = addcarry_u64(lo, t[i << 1], cc);
                let (d1, ee) = addcarry_u64(hi, t[(i << 1) + 1], ee);
                t[i << 1] = d0;
                t[(i << 1) + 1] = d1;
                cc = ee;
            }

            // Apply Montgomery reduction. We use the following facts:
            //  - upper half is necessarily less than p
            //  - set_montyred() accepts a full-limbs input and outputs a
            //    value of at most p
            //  - set_add() tolerates an input operand equal to p provided
            //    that the sum is less than 2*p
            self.0.copy_from_slice(&t[..N]);
            self.set_montyred();
            let mut y = Self([0u64; N]);
            y.0.copy_from_slice(&t[N..]);
            self.set_add(&y);
        }

        /// Compute the square of this value.
        pub fn square(self) -> Self {
            let mut r = self;
            r.set_square();
            r
        }

        /// Halve this value.
        pub fn set_half(&mut self) {
            let mm = (self.0[0] & 1).wrapping_neg();
            let mut cc = 0;
            for i in 0..(N - 1) {
                let (d, ee) = addcarry_u64(
                    (self.0[i] >> 1) | (self.0[i + 1] << 63),
                    mm & HALF_MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
            let (d, _) = addcarry_u64(self.0[N - 1] >> 1,
                mm & HALF_MODULUS[N - 1], cc);
            self.0[N - 1] = d;
        }

        /// Compute the half of this value.
        pub fn half(self) -> Self {
            let mut r = self;
            r.set_half();
            r
        }

        /// Double this value.
        pub fn set_mul2(&mut self) {
            // Double (as an integer) and subtract the modulus.
            let mut cc = 0;
            let mut tb = 0;
            for i in 0..N {
                let w = self.0[i];
                let t = (w << 1) | tb;
                tb = w >> 63;
                let (d, ee) = subborrow_u64(t, MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }

            // add back modulus if the result was negative
            let mm = tb.wrapping_sub(cc as u64);
            cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        /// Compute the sum of this value with itself.
        pub fn mul2(self) -> Self {
            let mut r = self;
            r.set_mul2();
            r
        }

        /// Triple this value.
        pub fn set_mul3(&mut self) {
            let r = self.mul2();
            *self += &r;
        }

        /// Compute the triple of this value.
        pub fn mul3(self) -> Self {
            let mut r = self;
            r.set_mul3();
            r
        }

        /// Quadruple this value.
        pub fn set_mul4(&mut self) {
            self.set_mul2();
            self.set_mul2();
        }

        /// Compute the quadruple of this value.
        pub fn mul4(self) -> Self {
            let mut r = self;
            r.set_mul4();
            r
        }

        /// Multiply this value by a small signed integer k.
        pub fn set_mul_small(&mut self, k: i32) {
            // Get the absolute value of the multiplier (but remember the sign).
            let sk = (k >> 31) as u32;
            let ak = ((k as u32) ^ sk).wrapping_sub(sk);

            // Do the product over integers.
            let (d, mut hi) = umull(self.0[0], ak as u64);
            self.0[0] = d;
            for i in 1..N {
                let (d, ee) = umull_add(self.0[i], ak as u64, hi);
                self.0[i] = d;
                hi = ee;
            }

            // We write:
            //    p = p1*2^m + p0   (modulus)
            //    x = x1*2^m + x0   (unreduced product)
            // with:
            //    2^31 <= p1 < 2^32
            //    0 <= p0 < 2^m
            //    0 <= x0 < 2^m
            // Since the current value x is the product of the input (less
            // than p) by a multiplier of at most 2^31, we know that:
            //    0 <= x < p*2^31 < 2^(63+m)
            //    0 <= x1 < 2^63.
            // We compute:
            //    b = floor(x1/p1)
            // Analysis shows that floor(x/p) = b, b-1 or b+1.
            //
            // We thus obtain b, then increment it (unless b == p1), then
            // subtract b*p from x; we then add back p repeatedly until a
            // non-negative result is obtained. At most two conditional
            // additions are needed to achieve that result.
            //
            // Division by p1 can be done with the Granlund-Montgomery method:
            //    https://dl.acm.org/doi/10.1145/773473.178249
            // (LLVM usually applies that method, but may fail to do so if for
            // instance optimizing for code size on some platforms, thus it is
            // best to apply the method explicitly so that constant-time code
            // is more reliably achieved.)

            // Extract top word of x.
            let bl = BITLEN & 63;
            let x1 = if bl == 0 {
                    (self.0[N - 1] >> 32) | (hi << 32)
                } else if bl < 32 {
                    (self.0[N - 1] << (32 - bl)) | (self.0[N - 2] >> (32 + bl))
                } else if bl == 32 {
                    self.0[N - 1]
                } else {
                    (hi << (96 - bl)) | (self.0[N - 1] >> (bl - 32))
                };

            // Compute b = floor(x1/p1).
            let (_, t) = umull(x1, P1DIV_M);
            let b = (x1.wrapping_sub(t) >> 1).wrapping_add(t) >> 31;

            // Add 1 to b, unless b == p1 (we cannot have b > p1).
            let b = b + (P1.wrapping_sub(b) >> 63);

            // Subtract b*p from x.
            let mut cc1 = 0;
            let mut cc2 = 0;
            for i in 0..N {
                let (d, ee) = umull_add(b, MODULUS[i], cc1);
                cc1 = ee;
                let (d, ee) = subborrow_u64(self.0[i], d, cc2);
                self.0[i] = d;
                cc2 = ee;
            }
            let (mut hi, _) = subborrow_u64(hi, cc1, cc2);

            // Add p (at most twice) as long as the value is negative.
            for _ in 0..2 {
                let m = sgnw(hi);
                let mut cc = 0;
                for i in 0..N {
                    let (d, ee) = addcarry_u64(self.0[i], m & MODULUS[i], cc);
                    self.0[i] = d;
                    cc = ee;
                }
                hi = hi.wrapping_add(cc as u64);
            }

            // We computed self*|k|; we must adjust for the sign of k.
            self.set_condneg(sk);
        }

        /// Compute the product of this value by a small (unsigned) integer k.
        pub fn mul_small(self, k: i32) -> Self {
            let mut r = self;
            r.set_mul_small(k);
            r
        }

        /// Set this value to either a or b, depending on whether the control
        /// word ctl is 0x00000000 or 0xFFFFFFFF, respectively.
        /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        pub fn set_select(&mut self, a: &Self, b: &Self, ctl: u32) {
            let c = (ctl as u64) | ((ctl as u64) << 32);
            for i in 0..N {
                let wa = a.0[i];
                let wb = b.0[i];
                self.0[i] = wa ^ (c & (wa ^ wb));
            }
        }

        /// Return a or b, if ctl is 0x00000000 or 0xFFFFFFFF, respectively.
        /// ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        pub fn select(a: &Self, b: &Self, ctl: u32) -> Self {
            let mut r = Self::ZERO;
            r.set_select(a, b, ctl);
            r
        }

        /// Set this value to rhs if ctl is 0xFFFFFFFF; leave it unchanged if
        /// ctl is 0x00000000.
        /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
            let c = (ctl as u64) | ((ctl as u64) << 32);
            for i in 0..N {
                let wa = self.0[i];
                let wb = rhs.0[i];
                self.0[i] = wa ^ (c & (wa ^ wb));
            }
        }

        /// Negate this value if ctl is 0xFFFFFFFF; leave it unchanged if
        /// ctl is 0x00000000.
        /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        pub fn set_condneg(&mut self, ctl: u32) {
            let v = -(self as &Self);
            self.set_cond(&v, ctl);
        }

        /// Exchange the values of a and b is ctl is 0xFFFFFFFF; leave both
        /// values unchanged if ctl is 0x00000000.
        /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
        pub fn condswap(a: &mut Self, b: &mut Self, ctl: u32) {
            let c = (ctl as u64) | ((ctl as u64) << 32);
            for i in 0..N {
                let wa = a.0[i];
                let wb = b.0[i];
                let wc = c & (wa ^ wb);
                a.0[i] = wa ^ wc;
                b.0[i] = wb ^ wc;
            }
        }

        // Set this value to (u*f+v*g)/2^64. Coefficients f
        // and g are provided as u64, but they are signed integers in the
        // [-2^62..+2^62] range.
        fn set_montylin(&mut self, u: &Self, v: &Self, f: u64, g: u64) {
            // Make sure f and g are non-negative.
            let sf = sgnw(f);
            let f = (f ^ sf).wrapping_sub(sf);
            let tu = Self::select(u, &-u, sf as u32); 
            let sg = sgnw(g);
            let g = (g ^ sg).wrapping_sub(sg);
            let tv = Self::select(v, &-v, sg as u32); 

            let (d, mut cc) = umull_x2(tu.0[0], f, tv.0[0], g);
            self.0[0] = d;
            for i in 1..N {
                let (d, hi) = umull_x2_add(tu.0[i], f, tv.0[i], g, cc);
                self.0[i] = d;
                cc = hi;
            }
            let up = cc;

            // Montgomery reduction (one round)
            let k = self.0[0].wrapping_mul(P0I);
            let (_, mut cc) = umull_add(k, MODULUS[0], self.0[0]);
            for i in 1..N {
                let (d, hi) = umull_add2(k, MODULUS[i], self.0[i], cc);
                self.0[i - 1] = d;
                cc = hi;
            }
            let (d, cc1) = addcarry_u64(up, cc, 0);
            self.0[N - 1] = d;

            // |f| <= 2^62 and |g| <= 2^62, therefore |u*f + v*g| < p*2^63
            // We added less than p*2^64, and divided by 2^64, so the result
            // is less than 2*p and a single conditional subtraction is enough.
            let mut cc2 = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(self.0[i], MODULUS[i], cc2);
                self.0[i] = d;
                cc2 = ee;
            }
            let mm = (cc1 as u64).wrapping_sub(cc2 as u64);
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = addcarry_u64(self.0[i], mm & MODULUS[i], cc);
                self.0[i] = d;
                cc = ee;
            }
        }

        fn montylin(a: &Self, b: &Self, f: u64, g: u64) -> Self {
            let mut r = Self::ZERO;
            r.set_montylin(a, b, f, g);
            r
        }

        // Set this value to abs((a*f + b*g)/2^31). Values a and b are
        // interpreted as plain integers (not modular). Coefficients f and
        // g are provided as u64 but they really are signed integers in the
        // [-2^31..+2^31] range (inclusive). The low 31 bits of a*f + b*g
        // are dropped (i.e. the division is assumed to be exact). The result
        // is assumed to fit in N limbs (extra high bits, if any, are
        // dropped). The absolute value of (a*f + b*g)/2^31 is computed.
        // Returned value is -1 (as a u64) if a*f + b*g was negative, 0
        // otherwise.
        fn set_lindiv31abs(&mut self, a: &Self, b: &Self, f: u64, g: u64)
            -> u64
        {
            // Replace f and g with abs(f) and abs(g), but remember the
            // original signs.
            let sf = sgnw(f);
            let f = (f ^ sf).wrapping_sub(sf);
            let sg = sgnw(g);
            let g = (g ^ sg).wrapping_sub(sg);

            // Compute a*f + b*g (upper word in 'up')
            let mut cc1 = 0;
            let mut cc2 = 0;
            let mut cc3 = 0;
            for i in 0..N {
                let (d1, ee1) = subborrow_u64(a.0[i] ^ sf, sf, cc1);
                cc1 = ee1;
                let (d2, ee2) = subborrow_u64(b.0[i] ^ sg, sg, cc2);
                cc2 = ee2;
                let (d3, hi3) = umull_x2_add(d1, f, d2, g, cc3);
                self.0[i] = d3;
                cc3 = hi3;
            }
            let up = cc3.wrapping_sub((cc1 as u64).wrapping_neg() & f)
                .wrapping_sub((cc2 as u64).wrapping_neg() & g);

            // Right-shift the result by 31 bits.
            for i in 0..(N - 1) {
                self.0[i] = (self.0[i] >> 31) | (self.0[i + 1] << 33);
            }
            self.0[N - 1] = (self.0[N - 1] >> 31) | (up << 33);

            // Negate the result if (a*f + b*g) was negative.
            let w = sgnw(up);
            let mut cc = 0;
            for i in 0..N {
                let (d, ee) = subborrow_u64(self.0[i] ^ w, w, cc);
                self.0[i] = d;
                cc = ee;
            }

            w
        }

        #[inline(always)]
        fn lindiv31abs(a: &Self, b: &Self, f: u64, g: u64) -> (Self, u64) {
            let mut r = Self::ZERO;
            let ng = r.set_lindiv31abs(a, b, f, g);
            (r, ng)
        }

        /// Divide this value by `y`. If `y` is zero, then this sets this value
        /// to zero.
        fn set_div(&mut self, y: &Self) {
            // a <- y
            // b <- p (modulus)
            // u <- x (self)
            // v <- 0
            //
            // Invariants:
            //    a*x = y*u mod p
            //    b*x = y*v mod p
            //    b is always odd
            //
            // At each step:
            //    if a is even, then:
            //        a <- a/2, u <- u/2 mod p
            //    else:
            //        if a < b:
            //            (a, u, b, v) <- (b, v, a, u)
            //        a <- (a - b)/2
            //        u <- (u - v)/2 mod p
            //
            // We optimize this algorithm following:
            //    https://eprint.iacr.org/2020/972

            let mut a = *y;
            let mut b = Self(MODULUS);
            let mut u = *self;
            let mut v = Self::ZERO;

            // Generic loop; each iteration reduces the sum of the sizes
            // of a and b by at least 31, and that sum starts at 2*BITLEN
            // (at most). We need to run it until the sum of the two lengths
            // is at most 64.
            const NUM1: usize = (2 * BITLEN - 34) / 31;
            for _ in 0..NUM1 {
                // Get approximations of a and b over 64 bits:
                //  - If len(a) <= 64 and len(b) <= 64, then we just
                //    use their values (low limbs).
                //  - Otherwise, with n = max(len(a), len(b)), we use:
                //       (a mod 2^31) + 2^31*floor(a / 2^(n - 33))
                //       (b mod 2^31) + 2^31*floor(b / 2^(n - 33))
                let mut c_hi = 0xFFFFFFFFFFFFFFFFu64;
                let mut c_lo = 0xFFFFFFFFFFFFFFFFu64;
                let mut a_hi = 0u64;
                let mut a_lo = 0u64;
                let mut b_hi = 0u64;
                let mut b_lo = 0u64;
                for j in (0..N).rev() {
                    let aw = a.0[j];
                    let bw = b.0[j];
                    a_hi ^= (a_hi ^ aw) & c_hi;
                    a_lo ^= (a_lo ^ aw) & c_lo;
                    b_hi ^= (b_hi ^ bw) & c_hi;
                    b_lo ^= (b_lo ^ bw) & c_lo;
                    c_lo = c_hi;
                    let mw = aw | bw;
                    c_hi &= ((mw | mw.wrapping_neg()) >> 63).wrapping_sub(1);
                }

                // If c_lo = 0, then we grabbed two words for a and b.
                // If c_lo != 0, then c_hi = 0 (they cannot be both non-zero
                // since that would mean that a = b = 0, but b is odd). In that
                // case, we grabbed one word (in a_hi and b_hi) and both values
                // fit in 64 bits.
                let s = lzcnt(a_hi | b_hi);
                let mut xa = (a_hi << s) | ((a_lo >> 1) >> (63 - s));
                let mut xb = (b_hi << s) | ((b_lo >> 1) >> (63 - s));
                xa = (xa & 0xFFFFFFFF80000000) | (a.0[0] & 0x000000007FFFFFFF);
                xb = (xb & 0xFFFFFFFF80000000) | (b.0[0] & 0x000000007FFFFFFF);

                // If c_lo != 0, then we should ignore the computed xa and xb,
                // and instead use the low limbs directly.
                xa ^= c_lo & (xa ^ a.0[0]);
                xb ^= c_lo & (xb ^ b.0[0]);

                // Compute the 31 inner iterations.
                let mut fg0 = 1u64;
                let mut fg1 = 1u64 << 32;
                for _ in 0..31 {
                    let a_odd = (xa & 1).wrapping_neg();
                    let (_, cc) = subborrow_u64(xa, xb, 0);
                    let swap = a_odd & (cc as u64).wrapping_neg();
                    let t1 = swap & (xa ^ xb);
                    xa ^= t1;
                    xb ^= t1;
                    let t2 = swap & (fg0 ^ fg1);
                    fg0 ^= t2;
                    fg1 ^= t2;
                    xa = xa.wrapping_sub(a_odd & xb);
                    fg0 = fg0.wrapping_sub(a_odd & fg1);
                    xa >>= 1;
                    fg1 <<= 1;
                }
                fg0 = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
                fg1 = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
                let f0 = (fg0 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g0 = (fg0 >> 32).wrapping_sub(0x7FFFFFFF);
                let f1 = (fg1 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g1 = (fg1 >> 32).wrapping_sub(0x7FFFFFFF);

                // Propagate updates to a, b, u and v.
                let (na, nega) = Self::lindiv31abs(&a, &b, f0, g0);
                let (nb, negb) = Self::lindiv31abs(&a, &b, f1, g1);
                let f0 = (f0 ^ nega).wrapping_sub(nega);
                let g0 = (g0 ^ nega).wrapping_sub(nega);
                let f1 = (f1 ^ negb).wrapping_sub(negb);
                let g1 = (g1 ^ negb).wrapping_sub(negb);
                let nu = Self::montylin(&u, &v, f0, g0);
                let nv = Self::montylin(&u, &v, f1, g1);
                a = na;
                b = nb;
                u = nu;
                v = nv;
            }

            const NUM2: usize = 2 * BITLEN - 31 * NUM1 - 2;

            // If y is non-zero, then the final GCD is 1, and
            // len(a) + len(b) <= NUM2 + 2 at this point (initially,
            // len(a) + len(b) <= 2*BITLEN, and each outer iteration reduces
            // the total by at least 31). Thus, the two values fit in one word
            // and we can finish the computation that way. We only need NUM2
            // iterations to reach the point where b = 1.
            let mut xa = a.0[0];
            let mut xb = b.0[0];
            let mut f0 = 1u64;
            let mut g0 = 0u64;
            let mut f1 = 0u64;
            let mut g1 = 1u64;
            for _ in 0..NUM2 {
                let a_odd = (xa & 1).wrapping_neg();
                let (_, cc) = subborrow_u64(xa, xb, 0);
                let swap = a_odd & (cc as u64).wrapping_neg();
                let t1 = swap & (xa ^ xb);
                xa ^= t1;
                xb ^= t1;
                let t2 = swap & (f0 ^ f1);
                f0 ^= t2;
                f1 ^= t2;
                let t3 = swap & (g0 ^ g1);
                g0 ^= t3;
                g1 ^= t3;
                xa = xa.wrapping_sub(a_odd & xb);
                f0 = f0.wrapping_sub(a_odd & f1);
                g0 = g0.wrapping_sub(a_odd & g1);
                xa >>= 1;
                f1 <<= 1;
                g1 <<= 1;
            }

            self.set_montylin(&u, &v, f1, g1);

            // If y != 0 then b = 1 at this point. If y == 0, then we
            // force the result to zero.
            let w = !y.iszero();
            let w = ((w as u64) << 32) | (w as u64);
            for i in 0..N {
                self.0[i] &= w;
            }

            // At this point, each outer iteration injected 31 extra doublings,
            // plus NUM2 for the last loop, for a total of NUM1*31 + NUM2.
            // Each montylin() call divided by 2^64, so in total we really
            // divided the value by 2^(64*(NUM1+1) - 31*NUM1 - NUM2).
            //
            // Moreover, both divisor and dividend were in Montgomery
            // representation, so the result is not in Montgomery representation
            // (the two R factors canceled each other). We want the result
            // in Montgomery representation, i.e. multiplied by 2^(64*N).
            // Therefore, we must multiply by 2^(33*NUM1 + 64 - NUM2 + 64*N),
            // which we need in
            self.set_mul(&Self::TFIXDIV);
        }

        pub fn set_invert(&mut self) {
            let r = *self;
            *self = Self::ONE;
            self.set_div(&r);
        }

        pub fn invert(self) -> Self {
            let mut r = Self::ONE;
            r.set_div(&self);
            r
        }

        /// Set this value to its square root. Returned value is 0xFFFFFFFF if
        /// the operation succeeded (value was indeed a quadratic residue), or
        /// 0x00000000 otherwise. On success, the chosen root is the one whose
        /// least significant bit (as an integer in [0..p-1]) is zero. On
        /// failure, this value is set to 0.
        pub fn set_sqrt(&mut self) -> u32 {
            const WIN_LEN: usize = 5;

            // Make a window.
            let mut ww = [*self; (1usize << WIN_LEN) - 1];
            for i in 1..ww.len() {
                if ((i + 1) & 1) == 0 {
                    ww[i] = ww[i >> 1].square();
                } else {
                    let z = &ww[i] * &ww[i - 1];
                    ww[i] = z;
                }
            }

            // Square and multiply algorithm, with exponent e = (p + 1)/4.
            // The exponent is not secret; we can do non-constant-time
            // lookups in the window, and omit multiplications for null digits.
            *self = ww[(SQRT_EH[SQRT_EH.len() - 1] as usize) - 1];
            for i in (0..(SQRT_EH.len() - 1)).rev() {
                for _ in 0..WIN_LEN {
                    self.set_square();
                }
                if SQRT_EH[i] != 0 {
                    self.set_mul(&ww[(SQRT_EH[i] as usize) - 1]);
                }
            }
            // Low 126 digits are all zero.
            for _ in 0..(WIN_LEN * SQRT_EL) {
                self.set_square();
            }

            // Check that the obtained value is indeed a square root of the
            // source value (which is still in ww[0]); if not, clear this
            // value.
            let r = self.square().equals(&ww[0]);
            let rw = (r as u64) | ((r as u64) << 32);
            for i in 0..N {
                self.0[i] &= rw;
            }

            // Conditionally negate this value, so that the chosen root
            // follows the expected convention.
            let ctl = ((self.encode()[0] as u32) & 1).wrapping_neg();
            self.set_condneg(ctl);

            r
        }

        /// Set this value to its fourth root. Returned value is 0xFFFFFFFF if
        /// the operation succeeded (value was indeed some element to the power of four), or
        /// 0x00000000 otherwise. On success, the chosen root is the one whose
        /// least significant bit (as an integer in [0..p-1]) is zero. On
        /// failure, this value is set to 0.
        pub fn set_fourth_root(&mut self) -> u32 {

            // TODO: make it to return the right root (as promised in the description)

            // TODO: move into a function and reuse in set_sqrt
            const WIN_LEN: usize = 5;
            // Make a window.
            let mut ww = [*self; (1usize << WIN_LEN) - 1];
            for i in 1..ww.len() {
                if ((i + 1) & 1) == 0 {
                    ww[i] = ww[i >> 1].square();
                } else {
                    let z = &ww[i] * &ww[i - 1];
                    ww[i] = z;
                }
            }

            // Square and multiply algorithm, with exponent e = (p + 1)/4.
            // The exponent is not secret; we can do non-constant-time
            // lookups in the window, and omit multiplications for null digits.
            *self = ww[(FOURTH_ROOT_EH[FOURTH_ROOT_EH.len() - 1] as usize) - 1];
            for i in (0..(FOURTH_ROOT_EH.len() - 1)).rev() {
                for _ in 0..WIN_LEN {
                    self.set_square();
                }
                if FOURTH_ROOT_EH[i] != 0 {
                    self.set_mul(&ww[(FOURTH_ROOT_EH[i] as usize) - 1]);
                }
            }
            // Low 126 digits are all zero.
            for _ in 0..(WIN_LEN * FOURTH_ROOT_EL) {
                self.set_square();
            }

            // Check that the obtained value is indeed a fourth root of the
            // source value (which is still in ww[0]); if not, clear this
            // value.
            let r = self.square().square().equals(&ww[0]);
            let rw = (r as u64) | ((r as u64) << 32);
            for i in 0..N {
                self.0[i] &= rw;
            }

            // Conditionally negate this value, so that the chosen root
            // follows the expected convention.
            let ctl = ((self.encode()[0] as u32) & 1).wrapping_neg();
            self.set_condneg(ctl);

            r
        }

        /// Compute the square root of this value. If this value is indeed a
        /// quadratic residue, then this returns (x, 0xFFFFFFFF), with x being
        /// the (unique) square root of this value whose least significant bit
        /// is zero (when normalized to an integer in [0..p-1]). If this value
        /// is not a quadratic residue, then this returns (zero, 0x00000000).
        pub fn sqrt(self) -> (Self, u32) {
            let mut x = self;
            let r = x.set_sqrt();
            (x, r)
        }

        /// Compute the fourth root of this value. If this value is indeed some
        /// element to the power of four, then this returns (x, 0xFFFFFFFF), with x being
        /// the (unique) fourth root of this value whose least significant bit
        /// is zero (when normalized to an integer in [0..p-1]). If this value
        /// is not some element to the power of four, then this returns (zero, 0x00000000).
        pub fn fourth_root(self) -> (Self, u32) {
            let mut x = self;
            let r = x.set_fourth_root();
            (x, r)
        }

        /// Raise this value to the power e. Exponent e is encoded in
        /// unsigned little-endian convention over exactly ebitlen bits.
        pub fn set_pow(&mut self, e: &[u8], ebitlen: usize) {
            self.set_pow_ext(e, 0, ebitlen);
        }

        /// Raise this value to the power e. Exponent e is encoded in
        /// unsigned little-endian convention, over exactly ebitlen bits,
        /// and starting at the bit offset eoff.
        pub fn set_pow_ext(&mut self, e: &[u8], eoff: usize, ebitlen: usize) {
            // TODO: implement a window optimization to make fewer
            // multiplications.
            let x = *self;
            *self = Self::ONE;
            for i in (eoff..(eoff + ebitlen)).rev() {
                let y = &*self * &x;
                let ctl = (((e[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                self.set_cond(&y, ctl);
                if i == eoff {
                    break;
                }
                self.set_square();
            }
        }

        /// Return this value to the power e (as a new element). Exponent e
        /// is encoded in unsigned little-endian convention over exactly
        /// ebitlen bits.
        pub fn pow(self, e: &[u8], ebitlen: usize) -> Self {
            let mut x = self;
            x.set_pow(e, ebitlen);
            x
        }

        /// Return this value to the power e (as a new element). Exponent e
        /// is encoded in unsigned little-endian convention over exactly
        /// ebitlen bits, and starting at the bit offset eoff.
        pub fn pow_ext(self, e: &[u8], eoff: usize, ebitlen: usize) -> Self {
            let mut x = self;
            x.set_pow_ext(e, eoff, ebitlen);
            x
        }

        /// Legendre symbol on this value. Return value is:
        ///   0   if this value is zero
        ///  +1   if this value is a non-zero quadratic residue
        ///  -1   if this value is not a quadratic residue
        pub fn legendre(self) -> i32 {
            // This is the same optimized binary GCD as in division, except
            // that we do not need to keep track of u and v. We can also
            // work directly on the Montgomery representation because R = 2^1184
            // is a square.
            let mut a = self;
            let mut b = Self(MODULUS);
            let mut ls = 0u64;

            // Generic loop; each iteration reduces the sum of the sizes
            // of a and b by at least 31, and that sum starts at 2*BITLEN
            // (at most). We need to run it until the sum of the two lengths
            // is at most 64.
            const NUM1: usize = (2 * BITLEN - 34) / 31;
            for _ in 0..NUM1 {
                // Get approximations of a and b over 64 bits:
                //  - If len(a) <= 64 and len(b) <= 64, then we just
                //    use their values (low limbs).
                //  - Otherwise, with n = max(len(a), len(b)), we use:
                //       (a mod 2^31) + 2^31*floor(a / 2^(n - 33))
                //       (b mod 2^31) + 2^31*floor(b / 2^(n - 33))
                let mut c_hi = 0xFFFFFFFFFFFFFFFFu64;
                let mut c_lo = 0xFFFFFFFFFFFFFFFFu64;
                let mut a_hi = 0u64;
                let mut a_lo = 0u64;
                let mut b_hi = 0u64;
                let mut b_lo = 0u64;
                for j in (0..N).rev() {
                    let aw = a.0[j];
                    let bw = b.0[j];
                    a_hi ^= (a_hi ^ aw) & c_hi;
                    a_lo ^= (a_lo ^ aw) & c_lo;
                    b_hi ^= (b_hi ^ bw) & c_hi;
                    b_lo ^= (b_lo ^ bw) & c_lo;
                    c_lo = c_hi;
                    let mw = aw | bw;
                    c_hi &= ((mw | mw.wrapping_neg()) >> 63).wrapping_sub(1);
                }

                // If c_lo = 0, then we grabbed two words for a and b.
                // If c_lo != 0, then c_hi = 0 (they cannot be both non-zero
                // since that would mean that a = b = 0, but b is odd). In that
                // case, we grabbed one word (in a_hi and b_hi) and both values
                // fit in 64 bits.
                let s = lzcnt(a_hi | b_hi);
                let mut xa = (a_hi << s) | ((a_lo >> 1) >> (63 - s));
                let mut xb = (b_hi << s) | ((b_lo >> 1) >> (63 - s));
                xa = (xa & 0xFFFFFFFF80000000) | (a.0[0] & 0x000000007FFFFFFF);
                xb = (xb & 0xFFFFFFFF80000000) | (b.0[0] & 0x000000007FFFFFFF);

                // If c_lo != 0, then we should ignore the computed xa and xb,
                // and instead use the low limbs directly.
                xa ^= c_lo & (xa ^ a.0[0]);
                xb ^= c_lo & (xb ^ b.0[0]);

                // First 29 inner iterations.
                let mut fg0 = 1u64;
                let mut fg1 = 1u64 << 32;
                for _ in 0..29 {
                    let a_odd = (xa & 1).wrapping_neg();
                    let (_, cc) = subborrow_u64(xa, xb, 0);
                    let swap = a_odd & (cc as u64).wrapping_neg();
                    ls ^= swap & ((xa & xb) >> 1);
                    let t1 = swap & (xa ^ xb);
                    xa ^= t1;
                    xb ^= t1;
                    let t2 = swap & (fg0 ^ fg1);
                    fg0 ^= t2;
                    fg1 ^= t2;
                    xa = xa.wrapping_sub(a_odd & xb);
                    fg0 = fg0.wrapping_sub(a_odd & fg1);
                    xa >>= 1;
                    fg1 <<= 1;
                    ls ^= xb.wrapping_add(2) >> 2;
                }

                // Compute the updated a and b (low words only) to get enough
                // bits for the next two iterations.
                let fg0z = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
                let fg1z = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
                let f0 = (fg0z & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g0 = (fg0z >> 32).wrapping_sub(0x7FFFFFFF);
                let f1 = (fg1z & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g1 = (fg1z >> 32).wrapping_sub(0x7FFFFFFF);
                let mut a0 = a.0[0].wrapping_mul(f0)
                    .wrapping_add(b.0[0].wrapping_mul(g0)) >> 29;
                let mut b0 = a.0[0].wrapping_mul(f1)
                    .wrapping_add(b.0[0].wrapping_mul(g1)) >> 29;
                for _ in 0..2 {
                    let a_odd = (xa & 1).wrapping_neg();
                    let (_, cc) = subborrow_u64(xa, xb, 0);
                    let swap = a_odd & (cc as u64).wrapping_neg();
                    ls ^= swap & ((a0 & b0) >> 1);
                    let t1 = swap & (xa ^ xb);
                    xa ^= t1;
                    xb ^= t1;
                    let t2 = swap & (fg0 ^ fg1);
                    fg0 ^= t2;
                    fg1 ^= t2;
                    let t3 = swap & (a0 ^ b0);
                    a0 ^= t3;
                    b0 ^= t3;
                    xa = xa.wrapping_sub(a_odd & xb);
                    fg0 = fg0.wrapping_sub(a_odd & fg1);
                    a0 = a0.wrapping_sub(a_odd & b0);
                    xa >>= 1;
                    fg1 <<= 1;
                    a0 >>= 1;
                    ls ^= b0.wrapping_add(2) >> 2;
                }

                // Propagate updates to a and b.
                fg0 = fg0.wrapping_add(0x7FFFFFFF7FFFFFFF);
                fg1 = fg1.wrapping_add(0x7FFFFFFF7FFFFFFF);
                let f0 = (fg0 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g0 = (fg0 >> 32).wrapping_sub(0x7FFFFFFF);
                let f1 = (fg1 & 0xFFFFFFFF).wrapping_sub(0x7FFFFFFF);
                let g1 = (fg1 >> 32).wrapping_sub(0x7FFFFFFF);

                // Propagate updates to a, b, u and v.
                let (na, nega) = Self::lindiv31abs(&a, &b, f0, g0);
                let (nb, _)    = Self::lindiv31abs(&a, &b, f1, g1);
                ls ^= nega & (nb.0[0] >> 1);
                a = na;
                b = nb;
            }

            const NUM2: usize = 2 * BITLEN - 31 * NUM1 - 2;

            // If y is non-zero, then the final GCD is 1, and
            // len(a) + len(b) <= NUM2 + 2 at this point (initially,
            // len(a) + len(b) <= 2*BITLEN, and each outer iteration reduces
            // the total by at least 31). Thus, the two values fit in one word
            // and we can finish the computation that way. We only need NUM2
            // iterations to reach the point where b = 1.
            let mut xa = a.0[0];
            let mut xb = b.0[0];
            for _ in 0..NUM2 {
                let a_odd = (xa & 1).wrapping_neg();
                let (_, cc) = subborrow_u64(xa, xb, 0);
                let swap = a_odd & (cc as u64).wrapping_neg();
                ls ^= swap & ((xa & xb) >> 1);
                let t1 = swap & (xa ^ xb);
                xa ^= t1;
                xb ^= t1;
                xa = xa.wrapping_sub(a_odd & xb);
                xa >>= 1;
                ls ^= xb.wrapping_add(2) >> 2;
            }

            // At this point, if the source value was not zero, then the low
            // bit of ls contains the QR status (0 = square, 1 = non-square),
            // which we need to convert to the expected value (+1 or -1).
            // If y == 0, then we return 0, per the API.
            let r = 1u32.wrapping_sub(((ls as u32) & 1) << 1);
            (r & !(self.iszero() as u32)) as i32
        }

        /// Encode this value into bytes. Encoding uses little-endian, has
        /// a fixed size (for a given field), and is canonical.
        pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
            let mut r = self;
            r.set_montyred();
            let mut d = [0u8; Self::ENCODED_LENGTH];
            for i in 0..(N - 1) {
                d[(i * 8)..(i * 8 + 8)].copy_from_slice(&r.0[i].to_le_bytes());
            }
            d[((N - 1) * 8)..].copy_from_slice(&(r.0[N - 1]
                .to_le_bytes()[..Self::ENCODED_LENGTH - (N - 1) * 8]));
            d
        }

        fn set_decode_nocheck(&mut self, buf: &[u8]) {
            for i in 0..(N - 1) {
                self.0[i] = u64::from_le_bytes(
                    *<&[u8; 8]>::try_from(&buf[(8 * i)..(8 * i + 8)]).unwrap());
            }
            let mut w = 0u64;
            for j in 0..(Self::ENCODED_LENGTH - (N - 1) * 8) {
                w |= (buf[(N - 1) * 8 + j] as u64) << (8 * j);
            }
            self.0[N - 1] = w;
        }

        /// Decode the provided bytes into a field element. Returned values
        /// are the element and 0xFFFFFFFF on success, or the zero element and
        /// 0x00000000 on failure. A failure is reported if the source slice
        /// does not have exactly the canonical encoding length of a field
        /// element (Self::ENCODED_LENGTH), or if the source encodes
        /// an integer which is not in the [0..(p-1)] range.
        pub fn decode(buf: &[u8]) -> (Self, u32) {
            if buf.len() != Self::ENCODED_LENGTH {
                return (Self::ZERO, 0);
            }

            // decode raw value
            let mut r = Self::ZERO;
            r.set_decode_nocheck(buf);

            // check that the source is canonical; clear if invalid
            let (_, mut cc) = subborrow_u64(r.0[0], MODULUS[0], 0);
            for i in 1..N {
                let (_, ee) = subborrow_u64(r.0[i], MODULUS[i], cc);
                cc = ee;
            }
            let m = (cc as u64).wrapping_neg();
            for i in 0..N {
                r.0[i] &= m;
            }

            // convert to Montgomery representation
            r.set_mul(&Self::R2);
            (r, m as u32)
        }

        /// Set this element by decoding the provided bytes. The source slice
        /// can have arbitrary length; the bytes are interpreted with the
        /// unsigned little-endian convention (no sign bit), and the resulting
        /// integer is reduced modulo the field modulus p. By definition, this
        /// function does not enforce canonicality of the source value.
        pub fn set_decode_reduce(&mut self, buf: &[u8]) {
            let mut n = buf.len();
            if n == 0 {
                *self = Self::ZERO;
                return;
            }

            let mut tmp = [0u8; Self::ENCODED_LENGTH];
            const CLEN: usize = 8 * (N - 1);
            let mut nn = n % CLEN;
            if nn == 0 {
                nn = CLEN;
            }
            n -= nn;
            tmp[..nn].copy_from_slice(&buf[n..]);
            self.set_decode_nocheck(&tmp);

            while n > 0 {
                n -= CLEN;
                tmp[..CLEN].copy_from_slice(&buf[n..(n + CLEN)]);
                let mut d = Self::ZERO;
                d.set_decode_nocheck(&tmp);
                self.set_mul(&Self::TDEC);
                self.set_add(&d);
            }

            self.set_mul(&Self::R2);
        }

        /// Decode the provided bytes into a field element. The source slice
        /// can have arbitrary length; the bytes are interpreted with the
        /// unsigned little-endian convention (no sign bit), and the resulting
        /// integer is reduced modulo the field modulus p. By definition, this
        /// function does not enforce canonicality of the source value.
        pub fn decode_reduce(buf: &[u8]) -> Self {
            let mut x = Self::ZERO;
            x.set_decode_reduce(buf);
            x
        }

        /// Set this structure to a random field element (indistinguishable
        /// from uniform generation).
        pub fn set_rand<T: CryptoRng + RngCore>(&mut self, rng: &mut T) {
            let mut tmp = [0u8; Self::ENCODED_LENGTH + 16];
            rng.fill_bytes(&mut tmp);
            self.set_decode_reduce(&tmp);
        }

        /// Return a new random field element (indistinguishable from
        /// uniform generation).
        pub fn rand<T: CryptoRng + RngCore>(rng: &mut T) -> Self {
            let mut x = Self::ZERO;
            x.set_rand(rng);
            x
        }

        /// Get the "hash" of the value (low 64 bits of the Montgomery
        /// representation).
        pub fn hashcode(self) -> u64 {
            self.0[0]
        }
    }

    // ========================================================================

    impl fmt::Display for Fp {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let v_bytes = self.encode();
            let v_big = BigInt::from_bytes_le(Sign::Plus, &v_bytes);
    
            write!(f, "{}", v_big.to_string())
        }
    }

    impl Add<Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn add(self, other: Fp) -> Fp {
            let mut r = self;
            r.set_add(&other);
            r
        }
    }

    impl Add<&Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn add(self, other: &Fp) -> Fp {
            let mut r = self;
            r.set_add(other);
            r
        }
    }

    impl Add<Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn add(self, other: Fp) -> Fp {
            let mut r = *self;
            r.set_add(&other);
            r
        }
    }

    impl Add<&Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn add(self, other: &Fp) -> Fp {
            let mut r = *self;
            r.set_add(other);
            r
        }
    }

    impl AddAssign<Fp> for Fp {
        #[inline(always)]
        fn add_assign(&mut self, other: Fp) {
            self.set_add(&other);
        }
    }

    impl AddAssign<&Fp> for Fp {
        #[inline(always)]
        fn add_assign(&mut self, other: &Fp) {
            self.set_add(other);
        }
    }

    impl Div<Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn div(self, other: Fp) -> Fp {
            let mut r = self;
            r.set_div(&other);
            r
        }
    }

    impl Div<&Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn div(self, other: &Fp) -> Fp {
            let mut r = self;
            r.set_div(other);
            r
        }
    }

    impl Div<Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn div(self, other: Fp) -> Fp {
            let mut r = *self;
            r.set_div(&other);
            r
        }
    }

    impl Div<&Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn div(self, other: &Fp) -> Fp {
            let mut r = *self;
            r.set_div(other);
            r
        }
    }

    impl DivAssign<Fp> for Fp {
        #[inline(always)]
        fn div_assign(&mut self, other: Fp) {
            self.set_div(&other);
        }
    }

    impl DivAssign<&Fp> for Fp {
        #[inline(always)]
        fn div_assign(&mut self, other: &Fp) {
            self.set_div(other);
        }
    }

    impl Mul<Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn mul(self, other: Fp) -> Fp {
            let mut r = self;
            r.set_mul(&other);
            r
        }
    }

    impl Mul<&Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn mul(self, other: &Fp) -> Fp {
            let mut r = self;
            r.set_mul(other);
            r
        }
    }

    impl Mul<Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn mul(self, other: Fp) -> Fp {
            let mut r = *self;
            r.set_mul(&other);
            r
        }
    }

    impl Mul<&Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn mul(self, other: &Fp) -> Fp {
            let mut r = *self;
            r.set_mul(other);
            r
        }
    }

    impl MulAssign<Fp> for Fp {
        #[inline(always)]
        fn mul_assign(&mut self, other: Fp) {
            self.set_mul(&other);
        }
    }

    impl MulAssign<&Fp> for Fp {
        #[inline(always)]
        fn mul_assign(&mut self, other: &Fp) {
            self.set_mul(other);
        }
    }

    impl Neg for Fp {
        type Output = Fp;

        #[inline(always)]
        fn neg(self) -> Fp {
            let mut r = self;
            r.set_neg();
            r
        }
    }

    impl Neg for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn neg(self) -> Fp {
            let mut r = *self;
            r.set_neg();
            r
        }
    }

    impl Sub<Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn sub(self, other: Fp) -> Fp {
            let mut r = self;
            r.set_sub(&other);
            r
        }
    }

    impl Sub<&Fp> for Fp {
        type Output = Fp;

        #[inline(always)]
        fn sub(self, other: &Fp) -> Fp {
            let mut r = self;
            r.set_sub(other);
            r
        }
    }

    impl Sub<Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn sub(self, other: Fp) -> Fp {
            let mut r = *self;
            r.set_sub(&other);
            r
        }
    }

    impl Sub<&Fp> for &Fp {
        type Output = Fp;

        #[inline(always)]
        fn sub(self, other: &Fp) -> Fp {
            let mut r = *self;
            r.set_sub(other);
            r
        }
    }

    impl SubAssign<Fp> for Fp {
        #[inline(always)]
        fn sub_assign(&mut self, other: Fp) {
            self.set_sub(&other);
        }
    }

    impl SubAssign<&Fp> for Fp {
        #[inline(always)]
        fn sub_assign(&mut self, other: &Fp) {
            self.set_sub(other);
        }
    }

    // ========================================================================

    /// GF(p^2) implementation.
    #[derive(Clone, Copy, Debug)]
    pub struct Fp2 {
        x0: Fp,
        x1: Fp,
    }

    impl Fp2 {

        pub const ZERO: Self = Self { x0: Fp::ZERO, x1: Fp::ZERO, };
        pub const ONE: Self = Self { x0: Fp::ONE, x1: Fp::ZERO, };
        pub const TWO: Self = Self { x0: Fp::TWO, x1: Fp::ZERO, };
        pub const THREE: Self = Self { x0: Fp::THREE, x1: Fp::ZERO, };
        pub const FOUR: Self = Self { x0: Fp::FOUR, x1: Fp::ZERO, };
        pub const MINUS_ONE: Self = Self { x0: Fp::MINUS_ONE, x1: Fp::ZERO, };

        pub const ENCODED_LENGTH: usize = 2 * Fp::ENCODED_LENGTH;

        /// Non-quadratic residue.
        pub const NQR: Self = Self { x0: Fp::NQR_RE, x1: Fp::ONE };

        pub fn new(re: &Fp, im: &Fp) -> Self {
            Self { x0: *re, x1: *im, }
        }

        pub fn iszero(self) -> u32 {
            self.x0.iszero() & self.x1.iszero()
        }

        pub fn equals(self, rhs: &Self) -> u32 {
            self.x0.equals(&rhs.x0) & self.x1.equals(&rhs.x1)
        }

        fn set_add(&mut self, rhs: &Self) {
            self.x0 += &rhs.x0;
            self.x1 += &rhs.x1;
        }

        fn set_sub(&mut self, rhs: &Self) {
            self.x0 -= &rhs.x0;
            self.x1 -= &rhs.x1;
        }

        pub fn set_neg(&mut self) {
            self.x0.set_neg();
            self.x1.set_neg();
        }

        pub fn set_conj(&mut self) {
            self.x1.set_neg();
        }

        pub fn conj(self) -> Self {
            Self { x0: self.x0, x1: -&self.x1, }
        }

        fn set_mul(&mut self, rhs: &Self) {
            // a <- x0*y0
            // b <- x1*y1
            // c <- (x0 + x1)*(y0 + y1)
            // (x0 + i*x1)*(y0 + i*y1) = (x0*y0 - x1*y1) + i*(x0*y1 + y0*x1)
            //                         = (a - b) + i*(c - a - b)
            let a = &self.x0 * &rhs.x0;
            let b = &self.x1 * &rhs.x1;
            let c = &(&self.x0 + &self.x1) * &(&rhs.x0 + &rhs.x1);
            self.x0 = a;
            self.x0 -= &b;
            self.x1 = c;
            self.x1 -= &a;
            self.x1 -= &b;
        }

        pub fn set_square(&mut self) {
            // (x0 + i*x1)^2 = (x0^2 - x1^2) + 2*i*(x0*x1)
            //               = (x0 + x1)*(x0 - x1) + i*(2*x0*x1)
            let a = &self.x0 + &self.x1;
            let b = &self.x0 - &self.x1;
            self.x1 *= &self.x0;
            self.x1.set_mul2();
            self.x0 = a;
            self.x0 *= &b;
        }

        pub fn square(self) -> Self {
            let mut r = self;
            r.set_square();
            r
        }

        pub fn set_half(&mut self) {
            self.x0.set_half();
            self.x1.set_half();
        }

        pub fn half(self) -> Self {
            let mut r = self;
            r.set_half();
            r
        }

        pub fn set_mul2(&mut self) {
            self.x0.set_mul2();
            self.x1.set_mul2();
        }

        pub fn mul2(self) -> Self {
            let mut r = self;
            r.set_mul2();
            r
        }

        pub fn set_mul3(&mut self) {
            self.x0.set_mul3();
            self.x1.set_mul3();
        }

        pub fn mul3(self) -> Self {
            let mut r = self;
            r.set_mul3();
            r
        }

        pub fn set_mul4(&mut self) {
            self.x0.set_mul4();
            self.x1.set_mul4();
        }

        pub fn mul4(self) -> Self {
            let mut r = self;
            r.set_mul4();
            r
        }

        pub fn set_mul_small(&mut self, k: i32) {
            self.x0.set_mul_small(k);
            self.x1.set_mul_small(k);
        }

        pub fn mul_small(self, k: i32) -> Self {
            let mut r = self;
            r.set_mul_small(k);
            r
        }

        pub fn set_select(&mut self, a: &Self, b: &Self, ctl: u32) {
            self.x0.set_select(&a.x0, &b.x0, ctl);
            self.x1.set_select(&a.x1, &b.x1, ctl);
        }

        pub fn select(a: &Self, b: &Self, ctl: u32) -> Self {
            Self {
                x0: Fp::select(&a.x0, &b.x0, ctl),
                x1: Fp::select(&a.x1, &b.x1, ctl),
            }
        }

        pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
            self.x0.set_cond(&rhs.x0, ctl);
            self.x1.set_cond(&rhs.x1, ctl);
        }

        pub fn set_condneg(&mut self, ctl: u32) {
            let y0 = -(&self.x0);
            let y1 = -(&self.x1);
            self.x0.set_cond(&y0, ctl);
            self.x1.set_cond(&y1, ctl);
        }

        pub fn condswap(a: &mut Self, b: &mut Self, ctl: u32) {
            Fp::condswap(&mut a.x0, &mut b.x0, ctl);
            Fp::condswap(&mut a.x1, &mut b.x1, ctl);
        }

        fn set_div(&mut self, rhs: &Self) {
            // 1/(x0 + i*x1) = (x0 - i*x1)/(x0^2 + x1^2)
            let mut z = rhs.x0.square();
            z += &rhs.x1.square();
            z.set_invert();
            let mut r = *rhs;
            r.x1.set_neg();
            r.x0 *= &z;
            r.x1 *= &z;
            self.set_mul(&r);
        }

        pub fn div(self, rhs: &Self) -> Self {
            let mut r = self;
            r.set_div(rhs);
            r
        }

        pub fn set_invert(&mut self) {
            // 1/(x0 + i*x1) = (x0 - i*x1)/(x0^2 + x1^2)
            let mut z = self.x0.square();
            z += &self.x1.square();
            z.set_invert();
            self.x0 *= &z;
            self.x1 *= &z;
            self.x1.set_neg();
        }

        pub fn invert(self) -> Self {
            let mut r = self;
            r.set_invert();
            r
        }

        /// Legendre symbol on this value. Return value is:
        ///   0   if this value is zero
        ///  +1   if this value is a non-zero quadratic residue
        ///  -1   if this value is not a quadratic residue
        pub fn legendre(self) -> i32 {
            // x = x0 + i*x1 is a square in GF(p^2) if and only if
            // x0^2 + x1^2 is a square in GF(p). Moreover, x0^2 + x1^2 is
            // zero if and only if x is zero.
            (self.x0.square() + self.x1.square()).legendre()
        }

        /// Set this value to its square root. Returned value is 0xFFFFFFFF if
        /// the operation succeeded (value was indeed a quadratic residue), or
        /// 0x00000000 otherwise. On success, the chosen root is the one whose
        /// sign is 0 (i.e. if the "real part" is non-zero, then it is an even
        /// integer; if the "real part" is zero, then the "imaginary part" is
        /// an even integer). On failure, this value is set to 0.
        pub fn set_sqrt(&mut self) -> u32 {
            // x^p = (x0 + i*x1)^p = x0 - i*x1  (Frobenius automorphism)
            // Thus: x^(p+1) = (x0 + i*x1)*(x0 - i*x1) = x0^2 + x1^2, which
            // is an element of GF(p). All elements of GF(p) are squares in
            // GF(p^2), but x0^2 + x1^2 is not necessarily a square in GF(p).
            //
            // Let conj(p) = x^p = x0 - i*x1. Note that conj() is analogous to
            // the conjugate in complex numbers. In particular:
            //    conj(a + b) = conj(a) + conj(b)
            //    conj(a * b) = conj(a) * conj(b)
            // This implies that conj(x) is a square if and only if x is a
            // square, and conj(sqrt(x)) = sqrt(conj(x)). Thus, if x is a
            // square, then:
            //    (sqrt(x)*conj(sqrt(x)))^2 = x*conj(x) = x0^2 + x1^2
            // But sqrt(x)*conj(sqrt(x)) is in GF(p); therefore, if x is a
            // square, then x0^2 + x1^2 must be a square in GF(p).
            //
            // Suppose that y = y0 + i*y1 such that y^2 = x. Then:
            //   y0^2 - y1^2 = x0
            //   2*y0*y1 = x1
            // If x1 = 0 then:
            //    if x0.legendre() >= 0 then y = sqrt(x0)
            //                          else y = i*sqrt(-x0)
            // else:
            //    y0 != 0 (necessarily) and y1 = x1 / (2*y0)
            //    Thus:
            //       y0^4 - x0*y0^2 - (x1^2)/4 = 0
            //    Discriminant is delta = x0^2 + x1^2, which is always a square
            //    (see above). Therefore:
            //       y0^2 = (x0 +/- sqrt(delta))/2
            //    We can thus compute (x0 + sqrt(delta))/2 and check its
            //    Legendre symbol; we subtract sqrt(delta) from it if it is
            //    not a square. We then extract y0 as a square root of the
            //    result, and compute y1 from it.
            //
            // Main cost is the two square roots in GF(p) (for delta and
            // for y0); Legendre symbols and inversions are vastly faster.

            // sqrt_delta <- sqrt(x0^2 + x1^2)
            let (sqrt_delta, r1) = (self.x0.square() + self.x1.square()).sqrt();
            // y0sq <- (x0 + sqrt(delta)) / 2
            let mut y0sq = (self.x0 + sqrt_delta).half();
            // If x1 = 0, then replace y0sq with x0
            let x1z = self.x1.iszero();
            y0sq.set_cond(&self.x0, x1z);
            // Get the Legendre symbol and set nqr to 0xFFFFFFFF when y0sq
            // is not a square
            let ls = y0sq.legendre();
            let nqr = (ls >> 1) as u32;
            // If not a square:
            //    if x1 = 0, then y0sq contains x0 and we want -x0
            //    if x1 != 0, then y0sq <- y0sq - sqrt(delta)
            y0sq.set_condneg(nqr & x1z);
            y0sq.set_cond(&(y0sq - sqrt_delta), nqr & !x1z);
            // Get the square root.
            let (mut y0, r2) = y0sq.sqrt();
            let r = r1 & r2;
            // Compute y1 = x1 / (2*y0).
            let mut y1 = self.x1 / y0.mul2();
            // If x1 = 0, then the square root worked, and y1 = 0 at this point;
            // we must still exchange y0 and y1 if x0 was not a square.
            Fp::condswap(&mut y0, &mut y1, nqr & x1z);
            // Result goes into this object. If there was a failure (r == 0),
            // then we must clear both x0 and x1.
            self.x0.set_select(&Fp::ZERO, &y0, r);
            self.x1.set_select(&Fp::ZERO, &y1, r);
            // Sign mangement: negate the result if needed.
            let x0odd = ((self.x0.encode()[0] as u32) & 1).wrapping_neg();
            let x1odd = ((self.x1.encode()[0] as u32) & 1).wrapping_neg();
            let x0z = self.x0.iszero();
            self.set_condneg(x0odd | (x0z & x1odd));
            r
        }

        pub fn sqrt(self) -> (Self, u32) {
            let mut y = self;
            let r = y.set_sqrt();
            (y, r)
        }

        /// Set this value to its fourth root. Returned value is 0xFFFFFFFF if
        /// the operation succeeded (value was indeed a quadratic residue), or
        /// 0x00000000 otherwise. On success, the chosen root is the one whose
        /// sign is 0 (i.e. if the "real part" is non-zero, then it is an even
        /// integer; if the "real part" is zero, then the "imaginary part" is
        /// an even integer). On failure, this value is set to 0.
        pub fn set_fourth_root(&mut self) -> u32 {
            // TODO: constant time
            
            let delta = self.x0.square() + self.x1.square();

            // TODO: error handling
            let mut n = delta.fourth_root().0;
            let disc = (n.square() + self.x0).half();

            // TODO: error handling
            let disc_sqrt = disc.sqrt().0;

            // TODO: we do not know which of n or -n is correct, test with legendre
            let mut y02 = (n + disc_sqrt).half();

            // is_square
            if !y02.legendre() >= 0 {
                y02 -= n;
                n = -n;
            }

            // TODO: error handling
            let y0 = y02.sqrt().0;

            // When we have (y02 + y02 - n) = 0
            // Then we have y0^2 = y1^2

            let gamma = y02 + y02 - n;

            if gamma.equals(&Fp::from_u32(0)) == 1 {
                self.x0 = y0;
                self.x1 = y0;
            }

            let y1 = self.x1 / (Fp::from_u32(4) * y0 * gamma);

            self.x0 = y0;
            self.x1 = y1;

            // TODO: fix return value
            return 1
        }

        pub fn fourth_root(self) -> (Self, u32) {
            let mut y = self;
            let r = y.set_fourth_root();
            (y, r)
        }

        /// Raise this value to the power e. Exponent e is encoded in
        /// unsigned little-endian convention over exactly ebitlen bits.
        pub fn set_pow(&mut self, e: &[u8], ebitlen: usize) {
            self.set_pow_ext(e, 0, ebitlen);
        }

        /// Raise this value to the power e. Exponent e is encoded in
        /// unsigned little-endian convention, over exactly ebitlen bits,
        /// and starting at the bit offset eoff.
        pub fn set_pow_ext(&mut self, e: &[u8], eoff: usize, ebitlen: usize) {
            // TODO: implement a window optimization to make fewer
            // multiplications.
            let x = *self;
            *self = Self::ONE;
            for i in (eoff..(eoff + ebitlen)).rev() {
                let y = &*self * &x;
                let ctl = (((e[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                self.set_cond(&y, ctl);
                if i == eoff {
                    break;
                }
                self.set_square();
            }
        }

        /// Return this value to the power e (as a new element). Exponent e
        /// is encoded in unsigned little-endian convention over exactly
        /// ebitlen bits.
        pub fn pow(self, e: &[u8], ebitlen: usize) -> Self {
            let mut x = self;
            x.set_pow(e, ebitlen);
            x
        }

        /// Return this value to the power e (as a new element). Exponent e
        /// is encoded in unsigned little-endian convention over exactly
        /// ebitlen bits, and starting at the bit offset eoff.
        pub fn pow_ext(self, e: &[u8], eoff: usize, ebitlen: usize) -> Self {
            let mut x = self;
            x.set_pow_ext(e, eoff, ebitlen);
            x
        }

        pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
            let mut r = [0u8; Self::ENCODED_LENGTH];
            r[..Fp::ENCODED_LENGTH].copy_from_slice(&self.x0.encode());
            r[Fp::ENCODED_LENGTH..].copy_from_slice(&self.x1.encode());
            r
        }

        pub fn decode(buf: &[u8]) -> (Self, u32) {
            if buf.len() != Self::ENCODED_LENGTH {
                return (Self::ZERO, 0);
            }
            let (mut x0, c0) = Fp::decode(&buf[..Fp::ENCODED_LENGTH]);
            let (mut x1, c1) = Fp::decode(&buf[Fp::ENCODED_LENGTH..]);
            let cx = c0 & c1;
            x0.set_cond(&Fp::ZERO, !cx);
            x1.set_cond(&Fp::ZERO, !cx);
            (Self { x0, x1, }, cx)
        }

        /// Set this structure to a random field element (indistinguishable
        /// from uniform generation).
        pub fn set_rand<T: CryptoRng + RngCore>(&mut self, rng: &mut T) {
            self.x0.set_rand(rng);
            self.x1.set_rand(rng);
        }

        /// Return a new random field element (indistinguishable from
        /// uniform generation).
        pub fn rand<T: CryptoRng + RngCore>(rng: &mut T) -> Self {
            let mut x = Self::ZERO;
            x.set_rand(rng);
            x
        }

        /// Set this structure to a random non-square field element
        /// (indistinguishable from uniform generation).
        pub fn set_rand_nonsquare<T: CryptoRng + RngCore>(&mut self,
            rng: &mut T)
        {
            // We get a random non-square by getting a random non-zero
            // value, squaring it, and multiplying by a known non-square.
            // The loop handles the case of getting zero randomly, which
            // is very unlikely in practice.
            loop {
                self.set_rand(rng);
                if self.iszero() != 0 {
                    continue;
                }
                break;
            }
            self.set_square();
            self.set_mul(&Self::NQR);
        }

        /// Return a new random non-square field element
        /// (indistinguishable from uniform generation).
        pub fn rand_nonsquare<T:CryptoRng + RngCore>(self, rng: &mut T)
            -> Self
        {
            let mut x = Self::ZERO;
            x.set_rand_nonsquare(rng);
            x
        }

        /// Raise this value to the power e. The exponent length (in bits)
        /// MUST be at most ebitlen. This is constant-time for both the
        /// base value (self) and the exponent (e); the exponent maximum
        /// size (ebitlen) is considered non-secret.
        pub fn set_pow_u64(&mut self, e: u64, ebitlen: usize) {
            match ebitlen {
                0 => {
                    *self = Self::ONE;
                }
                1 => {
                    self.set_cond(&Self::ONE, ((e as u32) & 1).wrapping_sub(1));
                }
                _ => {
                    let x = *self;
                    self.set_cond(&Self::ONE,
                        (((e >> (ebitlen - 1)) as u32) & 1).wrapping_sub(1));
                    for i in (0..(ebitlen - 1)).rev() {
                        self.set_square();
                        let y = &*self * &x;
                        self.set_cond(&y,
                            (((e >> i) as u32) & 1).wrapping_neg());
                    }
                }
            }
        }

        /// Return this value to the power e. The exponent length (in bits)
        /// MUST be at most ebitlen. This is constant-time for both the
        /// base value (self) and the exponent (e); the exponent maximum
        /// size (ebitlen) is considered non-secret.
        pub fn pow_u64(self, e: u64, ebitlen: usize) -> Self {
            let mut x = self;
            x.set_pow_u64(e, ebitlen);
            x
        }

        /// Raise this value to the power e. The exponent is considered
        /// non-secret.
        pub fn set_pow_u64_vartime(&mut self, e: u64) {
            match e {
                0 => {
                    *self = Self::ONE;
                }
                1 => {
                    return;
                }
                2 => {
                    self.set_square();
                }
                3 => {
                    *self *= self.square();
                }
                4 => {
                    self.set_square();
                    self.set_square();
                }
                _ => {
                    let xx = self.square();
                    let xw = [*self, xx, xx * &*self];
                    let mut j = 63 - e.leading_zeros();
                    j &= !1u32;
                    *self = xw[((e >> j) as usize) - 1];
                    while j > 0 {
                        j -= 2;
                        self.set_square();
                        self.set_square();
                        let k = ((e >> j) as usize) & 3;
                        if k > 0 {
                            self.set_mul(&xw[k - 1]);
                        }
                    }
                }
            }
        }

        /// Return this value to the power e. The exponent is considered
        /// non-secret.
        pub fn pow_u64_vartime(self, e: u64) -> Self {
            let mut x = self;
            x.set_pow_u64_vartime(e);
            x
        }

        /// Get the "hash" of the value. For x = x0 + i*x1, this is:
        ///    (hashcode(x0) << 1) | (hashcode(x1) & 1)
        /// i.e. bit 0 is bit 0 of x1, and bits 1..63 are bits 0..62 of x0
        /// (both in Montgomery representation).
        pub fn hashcode(self) -> u64 {
            (self.x0.hashcode() << 1) | (self.x1.hashcode() & 1)
        }
    }

    // ========================================================================

    impl fmt::Display for Fp2 {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let r = self.encode();

            let x0_bytes = &r[..Fp::ENCODED_LENGTH];
            let x1_bytes = &r[Fp::ENCODED_LENGTH..];

            let x0_big = BigInt::from_bytes_le(Sign::Plus, x0_bytes);
            let x1_big = BigInt::from_bytes_le(Sign::Plus, x1_bytes);

            write!(f, "i*{} + {}", x1_big, x0_big)
        }
    }

    impl Add<Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn add(self, other: Fp2) -> Fp2 {
            let mut r = self;
            r.set_add(&other);
            r
        }
    }

    impl Add<&Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn add(self, other: &Fp2) -> Fp2 {
            let mut r = self;
            r.set_add(other);
            r
        }
    }

    impl Add<Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn add(self, other: Fp2) -> Fp2 {
            let mut r = *self;
            r.set_add(&other);
            r
        }
    }

    impl Add<&Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn add(self, other: &Fp2) -> Fp2 {
            let mut r = *self;
            r.set_add(other);
            r
        }
    }

    impl AddAssign<Fp2> for Fp2 {
        #[inline(always)]
        fn add_assign(&mut self, other: Fp2) {
            self.set_add(&other);
        }
    }

    impl AddAssign<&Fp2> for Fp2 {
        #[inline(always)]
        fn add_assign(&mut self, other: &Fp2) {
            self.set_add(other);
        }
    }

    impl Div<Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn div(self, other: Fp2) -> Fp2 {
            let mut r = self;
            r.set_div(&other);
            r
        }
    }

    impl Div<&Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn div(self, other: &Fp2) -> Fp2 {
            let mut r = self;
            r.set_div(other);
            r
        }
    }

    impl Div<Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn div(self, other: Fp2) -> Fp2 {
            let mut r = *self;
            r.set_div(&other);
            r
        }
    }

    impl Div<&Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn div(self, other: &Fp2) -> Fp2 {
            let mut r = *self;
            r.set_div(other);
            r
        }
    }

    impl DivAssign<Fp2> for Fp2 {
        #[inline(always)]
        fn div_assign(&mut self, other: Fp2) {
            self.set_div(&other);
        }
    }

    impl DivAssign<&Fp2> for Fp2 {
        #[inline(always)]
        fn div_assign(&mut self, other: &Fp2) {
            self.set_div(other);
        }
    }

    impl Mul<Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn mul(self, other: Fp2) -> Fp2 {
            let mut r = self;
            r.set_mul(&other);
            r
        }
    }

    impl Mul<&Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn mul(self, other: &Fp2) -> Fp2 {
            let mut r = self;
            r.set_mul(other);
            r
        }
    }

    impl Mul<Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn mul(self, other: Fp2) -> Fp2 {
            let mut r = *self;
            r.set_mul(&other);
            r
        }
    }

    impl Mul<&Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn mul(self, other: &Fp2) -> Fp2 {
            let mut r = *self;
            r.set_mul(other);
            r
        }
    }

    impl MulAssign<Fp2> for Fp2 {
        #[inline(always)]
        fn mul_assign(&mut self, other: Fp2) {
            self.set_mul(&other);
        }
    }

    impl MulAssign<&Fp2> for Fp2 {
        #[inline(always)]
        fn mul_assign(&mut self, other: &Fp2) {
            self.set_mul(other);
        }
    }

    impl Neg for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn neg(self) -> Fp2 {
            let mut r = self;
            r.set_neg();
            r
        }
    }

    impl Neg for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn neg(self) -> Fp2 {
            let mut r = *self;
            r.set_neg();
            r
        }
    }

    impl Sub<Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn sub(self, other: Fp2) -> Fp2 {
            let mut r = self;
            r.set_sub(&other);
            r
        }
    }

    impl Sub<&Fp2> for Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn sub(self, other: &Fp2) -> Fp2 {
            let mut r = self;
            r.set_sub(other);
            r
        }
    }

    impl Sub<Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn sub(self, other: Fp2) -> Fp2 {
            let mut r = *self;
            r.set_sub(&other);
            r
        }
    }

    impl Sub<&Fp2> for &Fp2 {
        type Output = Fp2;

        #[inline(always)]
        fn sub(self, other: &Fp2) -> Fp2 {
            let mut r = *self;
            r.set_sub(other);
            r
        }
    }

    impl SubAssign<Fp2> for Fp2 {
        #[inline(always)]
        fn sub_assign(&mut self, other: Fp2) {
            self.set_sub(&other);
        }
    }

    impl SubAssign<&Fp2> for Fp2 {
        #[inline(always)]
        fn sub_assign(&mut self, other: &Fp2) {
            self.set_sub(other);
        }
    }

} } // End of macro: define_fp_core

pub(crate) use define_fp_core;

// ========================================================================

// Macro expectations:
#[cfg(test)]
macro_rules! define_fp_tests { () => {

    use super::{Fp, Fp2};
    use num_bigint::{BigInt, Sign, ToBigInt};
    use sha2::{Sha256, Digest};

    fn check_fp_ops(va: &[u8], vb: &[u8], with_sqrt: bool) {
        let mut zpww = [0u32; super::N * 2];
        for i in 0..super::N {
            zpww[2 * i] = super::MODULUS[i] as u32;
            zpww[2 * i + 1] = (super::MODULUS[i] >> 32) as u32;
        }
        let zp = BigInt::from_slice(Sign::Plus, &zpww);
        let zpz = &zp << 64;

        let a = Fp::decode_reduce(va);
        let b = Fp::decode_reduce(vb);
        let za = BigInt::from_bytes_le(Sign::Plus, va);
        let zb = BigInt::from_bytes_le(Sign::Plus, vb);

        let vc = a.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = &za % &zp;
        assert!(zc == zd);

        let c = a + b;
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za + &zb) % &zp;
        assert!(zc == zd);

        let c = a - b;
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = ((&zpz + &za) - (&zb % &zp)) % &zp;
        assert!(zc == zd);

        let c = -a;
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&zp - (&za % &zp)) % &zp;
        assert!(zc == zd);

        let c = a * b;
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za * &zb) % &zp;
        assert!(zc == zd);

        let c = a.square();
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za * &za) % &zp;
        assert!(zc == zd);

        let mut c = a.half();
        c += c;
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = &za % &zp;
        assert!(zc == zd);

        let c = a.mul2();
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za + &za) % &zp;
        assert!(zc == zd);

        let c = a.mul3();
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za + &za + &za) % &zp;
        assert!(zc == zd);

        let c = a.mul4();
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (&za + &za + &za + &za) % &zp;
        assert!(zc == zd);

        let k = u32::from_le_bytes(*<&[u8; 4]>::try_from(&vb[0..4]).unwrap()) as i32;
        let c = a.mul_small(k);
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = ((&za % &zp) * k + &zpz) % &zp;
        assert!(zc == zd);

        let c = Fp::from_i32(k);
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (k.to_bigint().unwrap() + &zp) % &zp;
        assert!(zc == zd);

        let k = k as u32;
        let c = Fp::from_u32(k);
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = k.to_bigint().unwrap();
        assert!(zc == zd);

        let k = u64::from_le_bytes(*<&[u8; 8]>::try_from(&vb[0..8]).unwrap()) as i64;
        let c = Fp::from_i64(k);
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = (k.to_bigint().unwrap() + &zp) % &zp;
        assert!(zc == zd);

        let k = k as u64;
        let c = Fp::from_u64(k);
        let vc = c.encode();
        let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
        let zd = k.to_bigint().unwrap();
        assert!(zc == zd);

        let c = a / b;
        if b.iszero() != 0 {
            assert!(c.iszero() != 0);
        } else {
            let c = c * b;
            let vc = c.encode();
            let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
            let zd = &za % &zp;
            assert!(zc == zd);
        }

        let c = a.square();
        if c.iszero() != 0 {
            assert!(c.legendre() == 0);
        } else {
            assert!(c.legendre() == 1);
            let c = -c;
            assert!(c.legendre() == -1);
        }

        if with_sqrt {
            let (c, r) = (a * a).sqrt();
            assert!(r == 0xFFFFFFFF);
            let vc = c.encode();
            let zc = BigInt::from_bytes_le(Sign::Plus, &vc);
            assert!(zc.bit(0) == false);
            let zc = (&zc * &zc) % &zp;
            let zd = (&za * &za) % &zp;
            assert!(zc == zd);
            if a.iszero() == 0 {
                let (c, r) = (-(a * a)).sqrt();
                assert!(c.iszero() == 0xFFFFFFFF);
                assert!(r == 0x00000000);
            }
        }
    }

    #[test]
    fn fp_ops() {
        let mut va = [0u8; (Fp::ENCODED_LENGTH + 64) & !31usize];
        let mut vb = [0u8; (Fp::ENCODED_LENGTH + 64) & !31usize];
        for i in 0..300 {
            let mut sh = Sha256::new();
            for j in 0..(va.len() >> 5) {
                sh.update(((256 * i + 8 * j + 0) as u64).to_le_bytes());
                va[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
            }
            for j in 0..(vb.len() >> 5) {
                sh.update(((256 * i + 8 * j + 1) as u64).to_le_bytes());
                vb[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
            }
            if i == 10 || i == 12 {
                va.fill(0);
            }
            if i == 11 || i == 12 {
                vb.fill(0);
            }
            check_fp_ops(&va, &vb, i < 30);
        }
    }

    fn check_fp2_ops(va: &[u8], vb: &[u8], with_sqrt: bool) {
        let mut zpww = [0u32; super::N * 2];
        for i in 0..super::N {
            zpww[2 * i] = super::MODULUS[i] as u32;
            zpww[2 * i + 1] = (super::MODULUS[i] >> 32) as u32;
        }
        let zp = BigInt::from_slice(Sign::Plus, &zpww);

        let alen = va.len() >> 1;
        let blen = vb.len() >> 1;

        let a0 = Fp::decode_reduce(&va[..alen]);
        let a1 = Fp::decode_reduce(&va[alen..]);
        let b0 = Fp::decode_reduce(&vb[..blen]);
        let b1 = Fp::decode_reduce(&vb[blen..]);
        let za0 = BigInt::from_bytes_le(Sign::Plus, &a0.encode());
        let za1 = BigInt::from_bytes_le(Sign::Plus, &a1.encode());
        let zb0 = BigInt::from_bytes_le(Sign::Plus, &b0.encode());
        let zb1 = BigInt::from_bytes_le(Sign::Plus, &b1.encode());
        let a = Fp2::new(&a0, &a1);
        let b = Fp2::new(&b0, &b1);

        let c = a + b;
        let vc = c.encode();
        let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
        let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
        let zd0 = (&za0 + &zb0) % &zp;
        let zd1 = (&za1 + &zb1) % &zp;
        assert!(zc0 == zd0 && zc1 == zd1);

        let c = a - b;
        let vc = c.encode();
        let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
        let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
        let zd0 = (&zp + &za0 - &zb0) % &zp;
        let zd1 = (&zp + &za1 - &zb1) % &zp;
        assert!(zc0 == zd0 && zc1 == zd1);

        let c = a * b;
        let vc = c.encode();
        let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
        let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
        let zd0 = (&zp + ((&za0 * &zb0) % &zp) - ((&za1 * &zb1) % &zp)) % &zp;
        let zd1 = ((&za0 * &zb1) + (&za1 * &zb0)) % &zp;
        assert!(zc0 == zd0 && zc1 == zd1);

        let c = a.square();
        let vc = c.encode();
        let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
        let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
        let zd0 = (&zp + ((&za0 * &za0) % &zp) - ((&za1 * &za1) % &zp)) % &zp;
        let zd1 = ((&za0 * &za1) + (&za1 * &za0)) % &zp;
        assert!(zc0 == zd0 && zc1 == zd1);

        let c = a / b;
        if b.iszero() != 0 {
            assert!(c.iszero() == 0xFFFFFFFF);
        } else {
            let c = c * b;
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            assert!(zc0 == za0 && zc1 == za1);
        }

        let c = b.invert();
        if b.iszero() != 0 {
            assert!(c.iszero() == 0xFFFFFFFF);
        } else {
            let c = c * b;
            assert!(c.equals(&Fp2::ONE) == 0xFFFFFFFF);
        }

        if with_sqrt {
            let e = a * a;
            let (c, r) = e.sqrt();
            assert!(r == 0xFFFFFFFF);
            assert!((c * c).equals(&e) == 0xFFFFFFFF);
            let vc = c.encode();
            let zc0 = BigInt::from_bytes_le(Sign::Plus, &vc[..Fp::ENCODED_LENGTH]);
            let zc1 = BigInt::from_bytes_le(Sign::Plus, &vc[Fp::ENCODED_LENGTH..]);
            assert!(zc0.bit(0) == false);
            if zc0.sign() == Sign::NoSign {
                assert!(zc1.bit(0) == false);
            }
            if a.iszero() == 0 {
                assert!(e.legendre() == 1);
                let e = a * a * Fp2::NQR;
                assert!(e.legendre() == -1);
                let (c, r) = e.sqrt();
                assert!(r == 0);
                assert!(c.iszero() == 0xFFFFFFFF);
            } else {
                assert!(e.legendre() == 0);
            }

            if a0.iszero() == 0 {
                let f = Fp2::new(&a0, &Fp::ZERO);
                let (c, r) = f.sqrt();
                assert!(r == 0xFFFFFFFF);
                assert!((c * c).equals(&f) == 0xFFFFFFFF);
                let g = -f;
                let (c, r) = g.sqrt();
                assert!(r == 0xFFFFFFFF);
                assert!((c * c).equals(&g) == 0xFFFFFFFF);
            }
        }
    }

    #[test]
    fn fp2_ops() {
        let mut va = [0u8; (2 * Fp::ENCODED_LENGTH + 64) & !31usize];
        let mut vb = [0u8; (2 * Fp::ENCODED_LENGTH + 64) & !31usize];
        for i in 0..100 {
            let mut sh = Sha256::new();
            for j in 0..(va.len() >> 5) {
                sh.update(((16 * i + 8 * j + 0) as u64).to_le_bytes());
                va[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
            }
            for j in 0..(vb.len() >> 5) {
                sh.update(((16 * i + 8 * j + 1) as u64).to_le_bytes());
                vb[(32 * j)..(32 * j + 32)].copy_from_slice(&sh.finalize_reset());
            }
            if i == 10 || i == 12 {
                va.fill(0);
            }
            if i == 11 || i == 12 {
                vb.fill(0);
            }
            check_fp2_ops(&va, &vb, i < 10);
        }
    }

} } // End of macro: define_fp_tests

#[cfg(test)]
pub(crate) use define_fp_tests;
