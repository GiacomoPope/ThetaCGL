use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// ========================================================================
// GF(p)

/// An element of GF(p).
#[derive(Clone, Copy, Debug)]
pub struct GFp(u64);

impl GFp {

    // IMPLEMENTATION NOTES:
    // ---------------------
    //
    // Let R = 2^64 mod p. Element x is represented by x*R mod p, in the
    // 0..p-1 range (Montgomery representation). Values are never outside
    // of that range.
    //
    // Everything is constant-time. There are specialized "Boolean"
    // functions such as iszero() and equals() that output a u64 which
    // happens, in practice, to have value 0 (for "false") or 2^64-1
    // (for "true").

    pub const ENCODED_LENGTH: usize = 8;

    pub const N: usize = 1;

    /// GF(p) modulus: p = 2^64 - 2^32 + 1
    /// GF(p) modulus: p = 2^64 - 257
    // pub const MOD: u64 = 0xFFFFFFFF00000001;
    // TODO: take C into account
    pub const MOD: u64 = 0xfffffffffffffeff;
    pub const MODULUS: [u64; GFp::N] = [
        0xFFFFFFFFFFFFFEFF
    ];

    /// Element 0 in GF(p).
    pub const ZERO: GFp = GFp::from_u64_reduce(0);

    /// Element 1 in GF(p).
    pub const ONE: GFp = GFp::from_u64_reduce(1);

    /// Element -1 in GF(p).
    pub const MINUS_ONE: GFp = GFp::from_u64_reduce(GFp::MOD - 1); 

    pub const R_VAL: u64 = 0x0000000000000101;
    pub const MINUS_R_VAL: u64 = 0xFFFFFFFFFFFFFDFE;
    const DR_VAL: u64 = 0x0000000000000202;
    pub const TR_VAL: u64 = 0x0000000000000303;
    pub const QR_VAL: u64 = 0x0000000000000404;
    pub const R2_VAL: u64 = 0x0000000000010201;
    pub const P0I: u64 = 18374966859414961921;
    pub const TFIXDIV_VAL: u64 = 0x0000000410181004;
    pub const TDEC_VAL: u64 = 0x0000000000000101;

    const WIN_LEN: i32 = 4;
    const SQRT_EL: i32 = 1;
    const SQRT_EH: [i32; 15] = [12, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 3];

    const FOURTH_ROOT_EL: i32 = 1;
    const FOURTH_ROOT_EH: [i32; 15] = [14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1];

    pub const TWO: Self = GFp::from_u64_reduce(GFp::DR_VAL);
    pub const THREE: Self = GFp::from_u64_reduce(GFp::TR_VAL);
    pub const FOUR: Self = GFp::from_u64_reduce(GFp::QR_VAL);
    pub const TDEC: Self = GFp::from_u64_reduce(GFp::TDEC_VAL);

    // 2^128 mod p.
    // const R2: u64 = 0xFFFFFFFE00000001;
    // TODO: take C into account, now computed for fixed value (258)
    const R2: u64 = 0x10201;
    const MU: u64 = 0xff00ff00ff00ff01;

    // p = 2^64 - C
    pub const C: u64 = 257;
    const P_PLUS_ONE_HALF: u64 = 0x7fffffffffffff80;

    // Montgomery reduction: given x <= p*2^64 - 1,
    // return x/2^64 mod p (in the 0 to p-1 range).
    #[inline(always)]
    const fn montyred(x: u128) -> u64 {
        let r0 = x as u64;
        let q = (GFp::MU as u128 * r0 as u128) as u64;

        let qp = q as u128 * GFp::MOD as u128;
        let (o, c1) = x.overflowing_add(qp);
        let mut r = o >> 64 as u64;
        r = r as u128 + (c1 as u128 * 1<<64 as u128);
        let r1 = (r % GFp::MOD as u128) as u64;

        r1
    } 

    /// Build a GF(p) element from a 64-bit integer. Returned values
    /// are (r, c). If the source value v is lower than the modulus,
    /// then r contains the value v as an element of GF(p), and c is
    /// equal to 0xFFFFFFFFFFFFFFFF; otherwise, r contains zero (in
    /// GF(p)) and c is 0.
    pub fn from_u64(v: u64) -> (GFp, u64) {
        // Computation of c is a constant-time lower-than operation:
        // If v < 2^63 then v < p and its high bit is 0.
        // If v >= 2^63 then its high bit is 1, and v < p if and only if
        // the high bit of z = v-p is 1 (since p >= 2^63 itself).
        let z = v.wrapping_sub(GFp::MOD);
        let c = ((v & !z) >> 63).wrapping_sub(1);
        (GFp::from_u64_reduce(v & c), c)
    }

    /// Build a GF(p) element from a 64-bit integer. The provided
    /// integer is implicitly reduced modulo p.
    #[inline(always)]
    pub const fn from_u64_reduce(v: u64) -> GFp {
        // R^2 = 2^64 - 2^33 + 1 mod p.
        // With v < 2^64, we have R**2 * v < 2^128 - 2^97 + 2^64, which is in
        // range of montyred().
        GFp(GFp::montyred((v as u128) * (GFp::R2 as u128)))
    }

    /// Get the element as an integer, normalized in the 0..p-1
    /// range.
    #[inline(always)]
    pub const fn to_u64(self) -> u64 {
        // Conversion back to normal representation is only a matter of
        // dividing by 2^64 modulo p, and that is exactly what montyred()
        // computes.
        GFp::montyred(self.0 as u128)
    }

    /// Addition in GF(p)
    #[inline(always)]
    const fn add(self, rhs: Self) -> Self {
        // We compute a + b = a - (p - b).
        let (x1, c1) = self.0.overflowing_sub(GFp::MOD - rhs.0);
        let t = c1 as u64 * GFp::C;
        GFp(x1.wrapping_sub(t as u64))
    }

    /// Subtraction in GF(p)
    #[inline(always)]
    const fn sub(self, rhs: Self) -> Self {
        // See montyred() for details on the subtraction.
        let (x1, c1) = self.0.overflowing_sub(rhs.0);
        let t = c1 as u64 * GFp::C;
        GFp(x1.wrapping_sub(t as u64))
    }

    /// Negation in GF(p)
    #[inline(always)]
    const fn neg(self) -> Self {
        GFp::ZERO.sub(self)
    }

    /// Halving in GF(p) (division by 2).
    #[inline(always)]
    pub const fn half(self) -> Self {
        // If x is even, then this returned x/2.
        // If x is odd, then this returns (x-1)/2 + (p+1)/2 = (x+p)/2.
        GFp((self.0 >> 1).wrapping_add(
            (self.0 & 1).wrapping_neg() & GFp::P_PLUS_ONE_HALF))
    }

    /// Doubling in GF(p) (multiplication by 2).
    #[inline(always)]
    pub const fn double(self) -> Self {
        self.add(self)
    }

    /// Multiplication in GF(p)
    #[inline(always)]
    const fn mul(self, rhs: Self) -> Self {
        // If x < p and y < p, then x*y <= (p-1)^2, and is thus in
        // range of montyred().
        GFp(GFp::montyred((self.0 as u128) * (rhs.0 as u128)))
    }

    /// Squaring in GF(p)
    #[inline(always)]
    pub const fn square(self) -> Self {
        self.mul(self)
    }

    /// Multiple squarings in GF(p): return x^(2^n)
    pub fn msquare(self, n: u32) -> Self {
        let mut x = self;
        for _ in 0..n {
            x = x.square();
        }
        x
    }

    /// Inversion in GF(p); if the input is zero, then zero is returned.
    pub fn invert(self) -> Self {
        // This uses Fermat's little theorem: 1/x = x^(p-2) mod p.
        // We have p-2 = 0xFFFFFFFEFFFFFFFF. In the instructions below,
        // we call 'xj' the value x^(2^j-1).
        // We call 'yj' the value x^(2^j).
        let x = self;
        let y1 = x.square(); // x^2
        let y2 = y1.square(); // x^4
        let y3 = y2.square(); // x^8

        let x2 = x * x.square(); // x^3
        let x4 = x2 * x2.msquare(2); // x^15
        let x5 = x * x4.square(); // x^31
        let x7 = x2 * x5.msquare(2); // x^(2^7 - 1)
        let x10 = x5 * x5.msquare(5); // x^(2^10 - 1)
        let x15 = x5 * x10.msquare(5); // x^(2^15 - 1)
        let x16 = x * x15.square(); // x^(2^16 - 1)
        let x31 = x15 * x16.msquare(15); // x^(2^31 - 1)
        let x32 = x * x31.square(); // x^(2^32 - 1)

        // 2^64 - 257 - 2 = (2^32 - 1) * 2^32 + 2^32 - 259 = (2^32 - 1) * 2^32 + (2^16 - 1) * 2^16 + 2^16 - 259
        // c1 = 2^8 - 3 = (2^4 - 1) * 2^4 + 13
        let c1 = x4.msquare(4) * y3 * y2 * x;
        // c2 = 2^16 - 259 = 2^16 - 2^8 - 3 = (2^7 - 1) * 2^9 + 2^8 - 3
        let c2 = x7.msquare(9) * c1;
        let c3 = x16.msquare(16);
        let c4 = x32.msquare(32);

        return c4 * c3 * c2;
    }

    fn div(self, rhs: Self) -> Self {
        self * rhs.invert()
    }

    /// Test of equality with zero; return value is 0xFFFFFFFFFFFFFFFF
    /// if this value is equal to zero, or 0 otherwise.
    #[inline(always)]
    pub const fn iszero(self) -> u64 {
        // Since values are always canonicalized internally, 0 in GF(p)
        // is always represented by the integer 0.
        // x == 0 if and only if both x and -x have their high bit equal to 0.
        !((((self.0 | self.0.wrapping_neg()) as i64) >> 63) as u64)
    }

    /// Test of equality with one; return value is 0xFFFFFFFFFFFFFFFF
    /// if this value is equal to one, or 0 otherwise.
    #[inline(always)]
    pub const fn isone(self) -> u64 {
        self.equals(GFp::ONE)
    }

    /// Test of equality with minus one; return value is 0xFFFFFFFFFFFFFFFF
    /// if this value is equal to -1 mod p, or 0 otherwise.
    #[inline(always)]
    pub const fn isminusone(self) -> u64 {
        self.equals(GFp::MINUS_ONE)
    }

    /// Test of equality between two GF(p) elements; return value is
    /// 0xFFFFFFFFFFFFFFFF if the two values are equal, or 0 otherwise.
    #[inline(always)]
    pub const fn equals(self, rhs: Self) -> u64 {
        // Since internal representation is canonical, we can simply
        // do a xor between the two operands, and then use the same
        // expression as iszero().
        let t = self.0 ^ rhs.0;
        !((((t | t.wrapping_neg()) as i64) >> 63) as u64)
    }

    /// Legendre symbol: return x^((p-1)/2) (as a GF(p) element).
    pub fn legendre_gfp(self) -> GFp {
        // (p-1)/2 = 0x7fffffffffffff7f
        let x = self;
        let y1 = x.square(); // x^2

        let x2 = x * x.square(); // x^3
        let x4 = x2 * x2.msquare(2); // x^15
        let x5 = x * x4.square(); // x^31
        let x8 = x4 * x4.msquare(4);
        let x10 = x5 * x5.msquare(5); // x^(2^10 - 1)
        let x15 = x5 * x10.msquare(5); // x^(2^15 - 1)
        let x16 = x * x15.square(); // x^(2^16 - 1)
        let x31 = x15 * x16.msquare(15); // x^(2^31 - 1)

        // 2**63 - 129 = (2**31 - 1) * 2**32 + 2**32 - 129
        // = (2**31 - 1) * 2**32 + (2**16 - 1) * 2**16 + 2**16 - 129
        // = (2**31 - 1) * 2**32 + (2**16 - 1) * 2**16 + (2**8 - 1) * 2**8 + 2**8 - 129
        // = (2**31 - 1) * 2**32 + (2**16 - 1) * 2**16 + (2**8 - 1) * 2**8 + 127

        // 2^7 - 1 = (2**5 - 1) * 2**2 + 3
        let c1 = x5.msquare(2) * y1 * x;
        let c2 = x8.msquare(8);
        let c3 = x16.msquare(16);
        let c4 = x31.msquare(32);

        c4 * c3 * c2 * c1
    }

    pub fn legendre(self) -> i32 {
        let l = self.legendre_gfp();
        let c1 = l.equals(GFp::ONE);
        let c2 = l.equals(GFp::MINUS_ONE);
        let cc1 = (c1 & 1) as i32;
        let cc2 = (c2 & 1) as i32;

        cc1 - cc2
    }

    pub fn sqrt(self) -> (Self, u64) { 
        // Make a window.
        let mut ww = [self; (1usize << GFp::WIN_LEN) - 1];
        for i in 1..ww.len() {
            if ((i + 1) & 1) == 0 {
                ww[i] = ww[i >> 1].square();
            } else {
                let z = ww[i] * ww[i - 1];
                ww[i] = z;
            }
        }

        // Square and multiply algorithm, with exponent e = (p + 1)/4.
        // The exponent is not secret; we can do non-constant-time
        // lookups in the window, and omit multiplications for null digits.
        let mut x = ww[(GFp::SQRT_EH[GFp::SQRT_EH.len() - 1] as usize) - 1];
        for i in (0..(GFp::SQRT_EH.len() - 1)).rev() {
            for _ in 0..GFp::WIN_LEN {
                x = x.square();
            }
            if GFp::SQRT_EH[i] != 0 {
                x = x.mul(ww[(GFp::SQRT_EH[i] as usize) - 1]);
            }
        }
        // Low 126 digits are all zero.
        for _ in 0..(GFp::WIN_LEN * GFp::SQRT_EL) {
            x = x.square()
        }

        // Check that the obtained value is indeed a square root of the
        // source value (which is still in ww[0]); if not, clear this
        // value.
        let r = x.square().equals(ww[0]);
        let rw = (r as u64) | ((r as u64) << 32);
        x.0 &= rw;

        // Conditionally negate this value, so that the chosen root
        // follows the expected convention.
        let ctl = ((self.encode()[0] as u64) & 1).wrapping_neg();
        x.set_condneg(ctl);

        (x, r)
    }

    /// Set this value to its fourth root. Returned value is 0xFFFFFFFF if
    /// the operation succeeded (value was indeed some element to the power of four), or
    /// 0x00000000 otherwise. On success, the chosen root is the one whose
    /// least significant bit (as an integer in [0..p-1]) is zero. On
    /// failure, this value is set to 0.
    /// TODO: return type
    pub fn fourth_root(self) -> (Self, u64) {
        // Make a window.
        let mut ww = [self; (1usize << GFp::WIN_LEN) - 1];
        for i in 1..ww.len() {
            if ((i + 1) & 1) == 0 {
                ww[i] = ww[i >> 1].square();
            } else {
                let z = ww[i] * ww[i - 1];
                ww[i] = z;
            }
        }

        // Square and multiply algorithm, with exponent e = (p + 1)/8.
        // The exponent is not secret; we can do non-constant-time
        // lookups in the window, and omit multiplications for null digits.
        let mut x = ww[(GFp::FOURTH_ROOT_EH[GFp::FOURTH_ROOT_EH.len() - 1] as usize) - 1];
        for i in (0..(GFp::FOURTH_ROOT_EH.len() - 1)).rev() {
            for _ in 0..GFp::WIN_LEN {
                x = x.square();
            }
            if GFp::FOURTH_ROOT_EH[i] != 0 {
                x = x.mul(ww[(GFp::FOURTH_ROOT_EH[i] as usize) - 1]);
            }
        }
        // Low 126 digits are all zero.
        for _ in 0..(GFp::WIN_LEN * GFp::FOURTH_ROOT_EL) {
            x = x.square();
        }

        // Check that the obtained value is indeed a fourth root of the
        // source value (which is still in ww[0]); if not, clear this
        // value.
        let r = x.square().square().equals(ww[0]);
        let rw = (r as u64) | ((r as u64) << 32);
        x.0 &= rw;

        // Conditionally negate this value, so that the chosen root
        // follows the expected convention.
        let ctl = ((self.encode()[0] as u64) & 1).wrapping_neg();
        x.set_condneg(ctl);

        (x, r)
    }

    /// Select a value: this function returns x0 if c == 0, or x1 if
    /// c == 0xFFFFFFFFFFFFFFFF.
    #[inline(always)]
    pub fn select(c: u64, x0: GFp, x1: GFp) -> GFp {
        GFp(x0.0 ^ (c & (x0.0 ^ x1.0)))
    }

    /// Negate this value.
    #[inline]
    pub fn set_neg(&mut self) {
        self.0 = self.neg().0;
    }

    /// Double this value.
    #[inline]
    pub fn set_mul2(&mut self) {
        let r = self.double();
        self.0 = r.0;
    }

    /// Compute the sum of this value with itself.
    #[inline(always)]
    pub fn mul2(self) -> Self {
        let mut r = self;
        r.set_mul2();
        r
    }

    /// Compute the quadruple of this value.
    #[inline(always)]
    pub fn mul4(self) -> Self {
        let mut r = self;
        r.set_mul4();
        r
    }

    /// Halve this value.
    #[inline]
    pub fn set_half(&mut self) {
        self.0 = self.half().0;
    }

    /// Quadruple this value.
    #[inline]
    pub fn set_mul4(&mut self) {
        self.set_mul2();
        self.set_mul2();
    }

    /// Multiply this value by 8
    #[inline]
    pub fn set_mul8(&mut self) {
        self.set_mul2();
        self.set_mul2();
        self.set_mul2();
    }

    /// Set this value to either a or b, depending on whether the control
    /// word ctl is 0x00000000 or 0xFFFFFFFF, respectively.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    /// TODO: ctl
    #[inline]
    pub fn set_select(&mut self, a: &Self, b: &Self, c: u64) {
        // let c = (ctl as u64) | ((ctl as u64) << 32);
        self.0 = GFp::select(c, *a, *b).0;
    }

    /// Set this value to rhs if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    /// TODO: fix desc (u32)
    #[inline]
    pub fn set_cond(&mut self, rhs: &Self, ctl: u64) {
        let wa = self.0;
        let wb = rhs.0;
        self.0 = wa ^ (ctl & (wa ^ wb));
    }

    /// Exchange the values of a and b is ctl is 0xFFFFFFFF; leave both
    /// values unchanged if ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    /// TODO: ctl
    #[inline]
    pub fn condswap(a: &mut Self, b: &mut Self, c: u64) {
        // let c = (ctl as u64) | ((ctl as u64) << 32);
        let wa = a.0;
        let wb = b.0;
        let wc = c & (wa ^ wb);
        a.0 = wa ^ wc;
        b.0 = wb ^ wc;
    }

    pub fn set_invert(&mut self) {
        self.0 = self.invert().0;
    }

    /// Negate this value if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    /// TODO: fix description (ctl: u32)
    #[inline]
    pub fn set_condneg(&mut self, ctl: u64) {
        // let c = (ctl as u64) | ((ctl as u64) << 32);
        let v = self.neg();
        self.set_cond(&v, ctl);
    }

    pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
        let r = self.to_u64();
        let mut d = [0u8; Self::ENCODED_LENGTH];
        d[0..].copy_from_slice(
            &(r.to_le_bytes()[..Self::ENCODED_LENGTH]),
        );
        d
    }

    pub fn decode(buf: &[u8]) -> (Self, u64) {
        if buf.len() != Self::ENCODED_LENGTH {
            return (GFp::ZERO, 0);
        }
        GFp::from_u64(u64::from_le_bytes(*<&[u8; 8]>::try_from(&buf[0.. 8]).unwrap()))
    }

    /// Get the "hash" of the value (64 bits of the Montgomery
    /// representation).
    pub fn hashcode(self) -> u64 {
        self.0
    }

    /// Decode the provided bytes. The source slice
    /// can have arbitrary length; the bytes are interpreted with the
    /// unsigned little-endian convention (no sign bit), and the resulting
    /// integer is reduced modulo the field modulus p. By definition, this
    /// function does not enforce canonicality of the source value.
    #[inline]
    pub fn decode_reduce(buf: &[u8]) -> Self {
        let mut n = buf.len();
        let mut x = Self::ZERO;
        if n == 0 {
            return x;
        }

        let mut tmp = [0u8; Self::ENCODED_LENGTH];
        const CLEN: usize = 8;
        let mut nn = n % CLEN;
        if nn == 0 {
            nn = CLEN;
        }
        n -= nn;
        tmp[..nn].copy_from_slice(&buf[n..]);
        (x, _) = GFp::decode(&tmp); // TODO

        while n > 0 {
            n -= CLEN;
            tmp[..CLEN].copy_from_slice(&buf[n..(n + CLEN)]);
            let (d, _) = GFp::decode(&tmp); // TODO
            x = x * Self::TDEC;
            x = x + d;
        }

        // x * GFp::from_u64_reduce(Self::R2)
        x
    }

}

// We implement all the needed traits to allow use of the arithmetic
// operators on GF(p) values.

impl Add<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn add(self, other: GFp) -> GFp {
        GFp::add(self, other)
    }
}

impl AddAssign<GFp> for GFp {
    #[inline(always)]
    fn add_assign(&mut self, other: GFp) {
        *self = GFp::add(*self, other);
    }
}

impl Sub<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn sub(self, other: GFp) -> GFp {
        GFp::sub(self, other)
    }
}

impl SubAssign<GFp> for GFp {
    #[inline(always)]
    fn sub_assign(&mut self, other: GFp) {
        *self = GFp::sub(*self, other);
    }
}

impl Neg for GFp {
    type Output = GFp;

    #[inline(always)]
    fn neg(self) -> GFp {
        GFp::neg(self)
    }
}

impl Mul<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn mul(self, other: GFp) -> GFp {
        GFp::mul(self, other)
    }
}

impl MulAssign<GFp> for GFp {
    #[inline(always)]
    fn mul_assign(&mut self, other: GFp) {
        *self = GFp::mul(*self, other);
    }
}

impl Div<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn div(self, other: GFp) -> GFp {
        GFp::div(self, other)
    }
}

impl DivAssign<GFp> for GFp {
    #[inline(always)]
    fn div_assign(&mut self, other: GFp) {
        *self = GFp::div(*self, other);
    }
}

// ========================================================================
// Unit tests.

#[cfg(test)]
mod tests {
    use super::GFp;

    // A custom PRNG; not cryptographically secure, but good enough
    // for tests.
    #[cfg(test)]
    struct PRNG(u128);

    #[cfg(test)]
    impl PRNG {
        // A: a randomly selected prime integer.
        // B: a randomly selected odd integer.
        const A: u128 = 87981536952642681582438141175044346919;
        const B: u128 = 331203846847999889118488772711684568729;

        // Get the next pseudo-random 64-bit integer.
        fn next_u64(&mut self) -> u64 {
            self.0 = PRNG::A.wrapping_mul(self.0).wrapping_add(PRNG::B);
            (self.0 >> 64) as u64
        }

        // Fill buf[] with pseudo-random bytes.
        fn next(&mut self, buf: &mut [u8]) {
            let mut acc: u64 = 0;
            for i in 0..buf.len() {
                if (i & 7) == 0 {
                    acc = self.next_u64();
                }
                buf[i] = acc as u8;
                acc >>= 8;
            }
        }
    }

    fn check_gfp_eq(a: GFp, r: u128) {
        assert!(a.to_u64() == (r % (GFp::MOD as u128)) as u64);
    }

    fn test_gfp_ops(a: u64, b: u64) {
        let x = GFp::from_u64_reduce(a);
        let y = GFp::from_u64_reduce(b);
        let wa = a as u128;
        let wb = b as u128;
        check_gfp_eq(x + y, wa + wb);
        check_gfp_eq(x - y, (wa + (GFp::MOD as u128) * 2) - wb);
        check_gfp_eq(-y, (GFp::MOD as u128) * 2 - wb);
        check_gfp_eq(x * y, wa * wb);
        check_gfp_eq(x.square(), wa * wa);
        if a == 0 || a == GFp::MOD {
            check_gfp_eq(x.invert(), 0);
        } else {
            check_gfp_eq(x * x.invert(), 1);
        }
        assert!(x.half().double().equals(x) == 0xFFFFFFFFFFFFFFFF);
    }

    #[test]
    fn gfp_ops() {
        for i in 0..10 {
            let v: u64 = (i as u64) + GFp::MOD - 5;
            if i <= 4 {
                let (x, c) = GFp::from_u64(v);
                assert!(c == 0xFFFFFFFFFFFFFFFF);
                assert!(x.to_u64() == v);
                let y = GFp::from_u64_reduce(v);
                assert!(y.to_u64() == v);
            } else {
                let v2 = v - GFp::MOD;
                let (x, c) = GFp::from_u64(v);
                assert!(c == 0);
                assert!(x.to_u64() == 0);
                let y = GFp::from_u64_reduce(v);
                assert!(y.to_u64() == v2);
            }
        }

        test_gfp_ops(0, 0);
        test_gfp_ops(0, 1);
        test_gfp_ops(1, 0);
        test_gfp_ops(1, 1);
        test_gfp_ops(0, 0xFFFFFFFFFFFFFFFF);
        test_gfp_ops(0xFFFFFFFFFFFFFFFF, 0);
        test_gfp_ops(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF);
        test_gfp_ops(0, 0xFFFFFFFF00000000);
        test_gfp_ops(0xFFFFFFFF00000000, 0);
        test_gfp_ops(0xFFFFFFFF00000000, 0xFFFFFFFF00000000);
        let mut prng = PRNG(0);
        for _ in 0..10000 {
            let a = prng.next_u64();
            let b = prng.next_u64();
            test_gfp_ops(a, b);
        }
        assert!(GFp::ZERO.legendre_gfp().iszero() == 0xFFFFFFFFFFFFFFFF);
        let (s0, c0) = GFp::ZERO.sqrt();
        check_gfp_eq(s0, 0);
        assert!(c0 == 0xFFFFFFFFFFFFFFFF);
        for _ in 0..1000 {
            let x = GFp::from_u64_reduce((prng.next_u64() >> 1) + 1).square();
            assert!(x.legendre_gfp().equals(GFp::ONE) == 0xFFFFFFFFFFFFFFFF);
            let (r1, c1) = x.sqrt();
            assert!(r1.square().equals(x) == 0xFFFFFFFFFFFFFFFF);
            assert!(c1 == 0xFFFFFFFFFFFFFFFF);
            let y = x * GFp::from_u64_reduce(7);
            assert!(y.legendre_gfp().equals(GFp::MINUS_ONE) == 0xFFFFFFFFFFFFFFFF);
            let (r2, c2) = y.sqrt();
            assert!(r2.iszero() == 0xFFFFFFFFFFFFFFFF);
            assert!(c2 == 0);
        }
    }

}