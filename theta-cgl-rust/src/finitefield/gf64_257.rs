use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand_core::{CryptoRng, RngCore};

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

    /// GF(p) modulus: p = 2^64 - 257
    pub const MOD: u64 = 0xfffffffffffffeff;
    pub const MODULUS: [u64; GFp::N] = [0xFFFFFFFFFFFFFEFF];

    /// Element 0 in GF(p).
    pub const ZERO: GFp = GFp::from_u64_reduce(0);

    /// Element 1 in GF(p).
    pub const ONE: GFp = GFp::from_u64_reduce(1);

    /// Element -1 in GF(p).
    pub const MINUS_ONE: GFp = GFp::from_u64_reduce(GFp::MOD - 1);

    pub const TWO: Self = GFp::from_u64_reduce(2);
    pub const THREE: Self = GFp::from_u64_reduce(3);
    pub const FOUR: Self = GFp::from_u64_reduce(4);

    const P_PLUS_ONE_HALF: u64 = 0x7fffffffffffff80;

    // TODO: experiment with other reductions
    // Reduction modulo p = 2^64 - 2^8 - 1
    #[inline(always)]
    const fn fp_reduction(x: u128) -> u64 {
        // x = x_lo + 2^64 * x_hi
        //   = x_lo + x_hi + 2^8*x_hi
        //
        // The resulting v will have 64 + 9 bits max
        let x_lo = (x as u64) as u128;
        let x_hi = x >> 64;
        let v: u128 = x_lo + x_hi + (x_hi << 8);

        // Now we need to fold in these pieces where cc < 2^9
        let v_lo = v as u64;
        let v_hi = (v >> 64) as u64;

        let (r, c) = v_lo.overflowing_sub(GFp::MOD - ((v_hi << 8) + v_hi));
        let adj = ((c as u16).wrapping_neg() & 257) as u64;
        r.wrapping_sub(adj)
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
        GFp(GFp::fp_reduction(v as u128))
    }

    /// Get the element as an integer, normalized in the 0..p-1
    /// range.
    #[inline(always)]
    pub const fn to_u64(self) -> u64 {
        // Conversion back to normal representation is only a matter of
        // dividing by 2^64 modulo p, and that is exactly what fp_reduction()
        // computes.
        GFp::fp_reduction(self.0 as u128)
    }

    /// Addition in GF(p)
    #[inline(always)]
    const fn add(self, rhs: Self) -> Self {
        // We compute a + b = a - (p - b).
        let (x1, c1) = self.0.overflowing_sub(GFp::MOD - rhs.0);
        let t = (c1 as u64) + ((c1 as u64) << 8);
        GFp(x1.wrapping_sub(t as u64))
    }

    /// Subtraction in GF(p)
    #[inline(always)]
    const fn sub(self, rhs: Self) -> Self {
        // See fp_reduction() for details on the subtraction.
        let (x1, c1) = self.0.overflowing_sub(rhs.0);
        let t = (c1 as u64) + ((c1 as u64) << 8);
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
        GFp((self.0 >> 1).wrapping_add((self.0 & 1).wrapping_neg() & GFp::P_PLUS_ONE_HALF))
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
        // range of fp_reduction().
        GFp(GFp::fp_reduction((self.0 as u128) * (rhs.0 as u128)))
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
        let x = self.exp_p_minus_three_div_four();
        self * x.msquare(2)
    }

    fn div(self, rhs: Self) -> Self {
        self * rhs.invert()
    }

    /// Test of equality with zero; return value is 0xFFFFFFFF
    /// if this value is equal to zero, or 0 otherwise.
    #[inline(always)]
    pub const fn iszero(self) -> u32 {
        // Since values are always canonicalized internally, 0 in GF(p)
        // is always represented by the integer 0.
        // x == 0 if and only if both x and -x have their high bit equal to 0.
        !((((self.0 | self.0.wrapping_neg()) as i64) >> 63) as u32)
    }

    /// Test of equality with one; return value is 0xFFFFFFFF
    /// if this value is equal to one, or 0 otherwise.
    #[inline(always)]
    pub const fn isone(self) -> u32 {
        self.equals(&GFp::ONE)
    }

    /// Test of equality with minus one; return value is 0xFFFFFFFF
    /// if this value is equal to -1 mod p, or 0 otherwise.
    #[inline(always)]
    pub const fn isminusone(self) -> u32 {
        self.equals(&GFp::MINUS_ONE)
    }

    /// Test of equality between two GF(p) elements; return value is
    /// 0xFFFFFFFF if the two values are equal, or 0 otherwise.
    #[inline(always)]
    pub const fn equals(self, rhs: &Self) -> u32 {
        // Since internal representation is canonical, we can simply
        // do a xor between the two operands, and then use the same
        // expression as iszero().
        let t = self.0 ^ rhs.0;
        !((((t | t.wrapping_neg()) as i64) >> 63) as u32)
    }

    // Compute a^((p-3)/4), used for inversion and Legendre
    // to compute a^(p-2) and a^((p-1)/2)
    #[inline(always)]
    fn exp_p_minus_three_div_four(self) -> GFp {
        let x = self;
        let x2 = x * x.square(); // x^(2^2 - 1)
        let x3 = x * x2.square(); // x^(2^3 - 1)
        let x6 = x3 * x3.msquare(3); // x^(2^6 - 1)
        let x12 = x6 * x6.msquare(6); // x^(2^12 - 1)
        let x24 = x12 * x12.msquare(12); // x^(2^24 - 1)
        let x48 = x24 * x24.msquare(24); // x^(2^48 - 1)
        let x54 = x6 * x48.msquare(6); // x^(2^54 - 1)
        let x55 = x * x54.square(); // x^(2^55 - 1)
        x6 * x55.msquare(7) // x^(2^62 - 2**7) + (2^6 - 1) =
                            // x^(2^62 - 2**6 - 1) = (p-3)/4
    }

    pub fn legendre(self) -> i32 {
        let mut l = self.exp_p_minus_three_div_four();
        l = self * l.square();

        let c1 = l.equals(&GFp::ONE);
        let c2 = l.equals(&GFp::MINUS_ONE);
        let cc1 = (c1 & 1) as i32;
        let cc2 = (c2 & 1) as i32;

        cc1 - cc2
    }

    pub fn sqrt(self) -> (Self, u32) {
        // We use that p = 3 mod 4 to compute the square root by
        // x^((p+1)/4) mod p. We have (p+1)/4 = 2^62 - 2^6
        // In the instructions below, we call 'xj' the value x^(2^j-1).
        let x = self;
        let x2 = x * x.square();
        let x4 = x2 * x2.msquare(2);
        let x6 = x2 * x4.msquare(2);
        let x7 = x * x6.square();
        let x14 = x7 * x7.msquare(7);
        let x28 = x14 * x14.msquare(14);
        let x56 = x28 * x28.msquare(28);
        let mut y = x56.msquare(6);

        let ctl = ((y.0 as u32) & 1).wrapping_neg();
        y.set_condneg(ctl);

        let r = y.square().equals(&self);
        let r_64 = (r as u64) | ((r as u64) << 32);
        y.0 &= r_64;

        (y, r)
    }

    /// Set this value to its fourth root. Returned value is 0xFFFFFFFF if
    /// the operation succeeded (value was indeed some element to the power of four), or
    /// 0x00000000 otherwise. On success, the chosen root is the one whose
    /// least significant bit (as an integer in [0..p-1]) is zero. On
    /// failure, this value is set to 0.
    pub fn fourth_root(self) -> (Self, u32) {
        let x = self;
        let x2 = x * x.square();
        let x4 = x2 * x2.msquare(2);
        let x6 = x2 * x4.msquare(2);
        let x7 = x * x6.square();
        let x14 = x7 * x7.msquare(7);
        let x28 = x14 * x14.msquare(14);
        let x56 = x28 * x28.msquare(28);
        let mut y = x56.msquare(5);

        // Check that the obtained value is indeed a fourth root of the
        // source value (which is still in ww[0]); if not, clear this
        // value.
        let r = y.square().square().equals(&x);
        let rw = (r as u64) | ((r as u64) << 32);
        y.0 &= rw;

        // Conditionally negate this value, so that the chosen root
        // follows the expected convention.
        let ctl = ((y.0 as u32) & 1).wrapping_neg();

        y.set_condneg(ctl);

        (y, r)
    }

    /// Select a value: this function returns x0 if c == 0, or x1 if
    /// c == 0xFFFFFFFF.
    #[inline(always)]
    pub fn select(x0: &GFp, x1: &GFp, c: u32) -> GFp {
        let c_64 = (c as u64) | ((c as u64) << 32);
        GFp(x0.0 ^ (c_64 & (x0.0 ^ x1.0)))
    }

    /// Negate this value.
    #[inline]
    pub fn set_neg(&mut self) {
        self.0 = self.neg().0;
    }

    /// Multiplication in GF(p) by a small integer (less than 2^31).
    #[inline(always)]
    pub fn mul_small(self, rhs: i32) -> Self {
        // Since the 'rhs' value is not in Montgomery representation,
        // we need to do a manual reduction instead.
        let x = (self.0 as u128) * (rhs as u128);
        let xl = x as u64;
        let xh = (x >> 64) as u64;

        // Since rhs <= 2^31 - 1, we have xh <= 2^31 - 2, and
        // p - xh >= 2^64 - 2^31 - 2^8 - 3, which is close to 2^64;
        // thus, even if xl was not lower than p, the subtraction
        // will bring back the value in the proper range, and the
        // normal subtraction in GF(p) yields the proper result.
        let (r, c) = xl.overflowing_sub(GFp::MOD - ((xh << 8) + xh));
        let adj = ((c as u16).wrapping_neg() & 257) as u64;
        GFp(r.wrapping_sub(adj))
    }

    pub fn set_mul_small(&mut self, rhs: i32) {
        let r = self.mul_small(rhs);
        self.0 = r.0;
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
    #[inline]
    pub fn set_select(&mut self, a: &Self, b: &Self, ctl: u32) {
        self.0 = GFp::select(a, b, ctl).0;
    }

    /// Set this value to rhs if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline]
    pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
        let wa = self.0;
        let wb = rhs.0;
        let ctl_64 = (ctl as u64) | ((ctl as u64) << 32);
        self.0 = wa ^ (ctl_64 & (wa ^ wb));
    }

    /// Exchange the values of a and b is ctl is 0xFFFFFFFF; leave both
    /// values unchanged if ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline]
    pub fn condswap(a: &mut Self, b: &mut Self, c: u32) {
        let wa = a.0;
        let wb = b.0;
        let c_64 = (c as u64) | ((c as u64) << 32);
        let wc = c_64 & (wa ^ wb);
        a.0 = wa ^ wc;
        b.0 = wb ^ wc;
    }

    pub fn set_invert(&mut self) {
        self.0 = self.invert().0;
    }

    /// Negate this value if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline]
    pub fn set_condneg(&mut self, ctl: u32) {
        let v = self.neg();
        self.set_cond(&v, ctl);
    }

    pub fn encode(self) -> [u8; Self::ENCODED_LENGTH] {
        let r = self.to_u64();
        let mut d = [0u8; Self::ENCODED_LENGTH];
        d[0..].copy_from_slice(&(r.to_le_bytes()[..Self::ENCODED_LENGTH]));
        d
    }

    pub fn decode(buf: &[u8]) -> (Self, u64) {
        if buf.len() != Self::ENCODED_LENGTH {
            return (GFp::ZERO, 0);
        }
        GFp::from_u64(u64::from_le_bytes(
            *<&[u8; 8]>::try_from(&buf[0..8]).unwrap(),
        ))
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
        (x, _) = GFp::decode(&tmp);

        while n > 0 {
            n -= CLEN;
            tmp[..CLEN].copy_from_slice(&buf[n..(n + CLEN)]);
            let (d, _) = GFp::decode(&tmp);
            x = x + d;
        }

        x
    }

    pub fn set_decode_reduce(&mut self, buf: &[u8]) {
        let r = Self::decode_reduce(buf);
        self.0 = r.0;
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

impl Add<&GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn add(self, other: &GFp) -> GFp {
        GFp::add(self, *other)
    }
}

impl Add<GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn add(self, other: GFp) -> GFp {
        GFp::add(*self, other)
    }
}

impl Add<&GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn add(self, other: &GFp) -> GFp {
        GFp::add(*self, *other)
    }
}

impl AddAssign<GFp> for GFp {
    #[inline(always)]
    fn add_assign(&mut self, other: GFp) {
        *self = GFp::add(*self, other);
    }
}

impl AddAssign<&GFp> for GFp {
    #[inline(always)]
    fn add_assign(&mut self, other: &GFp) {
        *self = GFp::add(*self, *other);
    }
}

impl Sub<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn sub(self, other: GFp) -> GFp {
        GFp::sub(self, other)
    }
}

impl Sub<&GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn sub(self, other: &GFp) -> GFp {
        GFp::sub(self, *other)
    }
}

impl Sub<GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn sub(self, other: GFp) -> GFp {
        GFp::sub(*self, other)
    }
}

impl Sub<&GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn sub(self, other: &GFp) -> GFp {
        GFp::sub(*self, *other)
    }
}

impl SubAssign<GFp> for GFp {
    #[inline(always)]
    fn sub_assign(&mut self, other: GFp) {
        *self = GFp::sub(*self, other);
    }
}

impl SubAssign<&GFp> for GFp {
    #[inline(always)]
    fn sub_assign(&mut self, other: &GFp) {
        *self = GFp::sub(*self, *other);
    }
}

impl Neg for GFp {
    type Output = GFp;

    #[inline(always)]
    fn neg(self) -> GFp {
        GFp::neg(self)
    }
}

impl Neg for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn neg(self) -> GFp {
        GFp::neg(*self)
    }
}

impl Mul<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn mul(self, other: GFp) -> GFp {
        GFp::mul(self, other)
    }
}

impl Mul<&GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn mul(self, other: &GFp) -> GFp {
        GFp::mul(self, *other)
    }
}

impl Mul<GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn mul(self, other: GFp) -> GFp {
        GFp::mul(*self, other)
    }
}

impl Mul<&GFp> for &GFp {
    type Output = GFp;

    #[inline(always)]
    fn mul(self, other: &GFp) -> GFp {
        GFp::mul(*self, *other)
    }
}

impl MulAssign<GFp> for GFp {
    #[inline(always)]
    fn mul_assign(&mut self, other: GFp) {
        *self = GFp::mul(*self, other);
    }
}

impl MulAssign<&GFp> for GFp {
    #[inline(always)]
    fn mul_assign(&mut self, other: &GFp) {
        *self = GFp::mul(*self, *other);
    }
}

impl Div<GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn div(self, other: GFp) -> GFp {
        GFp::div(self, other)
    }
}

impl Div<&GFp> for GFp {
    type Output = GFp;

    #[inline(always)]
    fn div(self, other: &GFp) -> GFp {
        GFp::div(self, *other)
    }
}

impl DivAssign<GFp> for GFp {
    #[inline(always)]
    fn div_assign(&mut self, other: GFp) {
        *self = GFp::div(*self, other);
    }
}

impl DivAssign<&GFp> for GFp {
    #[inline(always)]
    fn div_assign(&mut self, other: &GFp) {
        *self = GFp::div(*self, *other);
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
        assert!(x.half().double().equals(&x) == 0xFFFFFFFF);
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
        assert!(GFp::ZERO.legendre() == 0);
        let (s0, c0) = GFp::ZERO.sqrt();
        check_gfp_eq(s0, 0);
        assert!(c0 == 0xFFFFFFFF);
        for _ in 0..1000 {
            let x = GFp::from_u64_reduce((prng.next_u64() >> 1) + 1).square();
            assert!(x.legendre() == 1);
            let (r1, c1) = x.sqrt();
            assert!(r1.square().equals(&x) == 0xFFFFFFFF);
            assert!(c1 == 0xFFFFFFFF);
            let y = x * GFp::from_u64_reduce(7);
            assert!(y.legendre() == -1);
            let (r2, c2) = y.sqrt();
            assert!(r2.iszero() == 0xFFFFFFFF);
            assert!(c2 == 0);

            let z = x.square();
            let (r3, c3) = z.fourth_root();
            assert!(r3.msquare(2).equals(&z) == 0xFFFFFFFF);
            assert!(c3 == 0xFFFFFFFF);
        }
    }
}
