use crate::finitefield::utils64::umull_add;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
// ========================================================================
// GF(p)

/// An element of GF(p).
#[derive(Clone, Copy, Debug)]
pub struct GF257(u64);

impl GF257 {
    // IMPLEMENTATION NOTES:
    // ---------------------
    //
    // Let R = 2^64 mod p. Element x is represented by x*R mod p, in the
    // 0..p-1 range (Montgomery representation). Values are never outside
    // of that range.
    //
    // Everything is constant-time. There are specialized "Boolean"
    // functions such as iszero() and equals() that output a u32 which
    // happens, in practice, to have value 0 (for "false") or 2^32-1
    // (for "true").

    /// GF(p) modulus: p = 2^64 - 2^8 - 1
    pub const MOD: u64 = 0xFFFFFFFFFFFFFEFF;

    /// Element 0 in GF(p).
    pub const ZERO: GF257 = GF257::from_u64_reduce(0);

    /// Element 1 in GF(p).
    pub const ONE: GF257 = GF257::from_u64_reduce(1);

    /// Element -1 in GF(p).
    pub const MINUS_ONE: GF257 = GF257::from_u64_reduce(GF257::MOD - 1);

    // 2^128 mod p.
    const R2: u64 = 0x0000000000010201;
    const MU: u64 = 0xff00ff00ff00ff01;

    // TODO: these are extra constants I needed for the extension
    // construction... Maybe these should be shifted to a specialised
    // extension...
    pub const ENCODED_LENGTH: usize = 8;
    pub const N: usize = 1; // Number of limbs
    pub const TWO: GF257 = GF257::from_u64_reduce(2);
    pub const THREE: GF257 = GF257::from_u64_reduce(3);
    pub const FOUR: GF257 = GF257::from_u64_reduce(4);
    pub const THREE_INV: GF257 = GF257::from_u64_reduce(0x5555555555555500);

    // Montgomery reduction: given x <= p*2^64 - 1,
    // return x/2^64 mod p (in the 0 to p-1 range).
    // #[inline(always)]
    // const fn montyred(x: u128) -> u64 {
    //     let r0 = x as u64;
    //     let q = r0.wrapping_mul(GF257::MU);

    //     let qp = q as u128 * GF257::MOD as u128;
    //     let (o, c1) = x.overflowing_add(qp);
    //     let mut r = o >> 64 as u64;
    //     r = r as u128 + ((c1 as u128) << 64 as u128);
    //     let r1 = (r % GF257::MOD as u128) as u64;

    //     r1
    // }

    #[inline(always)]
    const fn montyred(x: u128) -> u64 {
        // Split x into low and hi parts
        let x0 = x as u64;
        let x1 = (x >> 64) as u64;

        // Compute (mu * x0) mod 2^64 where
        // mu = -p mod 2^64
        let q = x0.wrapping_mul(GF257::MU);

        // We want to compute t where q*p + x = 2^64*t
        // So we only need the high 64 bits of x0 + q*p
        // and then we can add to this result x1
        let (_, t) = umull_add(q, GF257::MOD, x0);
        let (t, c) = t.overflowing_add(x1);

        // If the above addition overflowed we subtract p by
        // adding 2^8 + 1
        let out = t + ((c as u64) << 8) + (c as u64);

        // We're now in the case where out is smaller than 2^64
        // but we may still have out >= p. To ensure output is in
        // [0, p-1] we need another conditional subtraction
        // TODO: should we allow working with intermediate values
        // [0, 2^64 - 1]?

        // TODO: I don't like this >= check, should be done with sub.
        let c = out >= GF257::MOD;
        out.wrapping_add(((c as u64) << 8) + (c as u64))
    }

    /// Build a GF(p) element from a 64-bit integer. Returned values
    /// are (r, c). If the source value v is lower than the modulus,
    /// then r contains the value v as an element of GF(p), and c is
    /// equal to 0xFFFFFFFF; otherwise, r contains zero (in
    /// GF(p)) and c is 0.
    pub fn from_u64(v: u64) -> (GF257, u32) {
        // Computation of c is a constant-time lower-than operation:
        // If v < 2^63 then v < p and its high bit is 0.
        // If v >= 2^63 then its high bit is 1, and v < p if and only
        // the high bit of z = v-p is 1 (since p >= 2^63 itself).
        let z = v.wrapping_sub(GF257::MOD);
        let c = ((v & !z) >> 63).wrapping_sub(1);
        (GF257::from_u64_reduce(v & c), c as u32)
    }

    /// Build a GF(p) element from a 64-bit integer. The provided
    /// integer is implicitly reduced modulo p.
    #[inline(always)]
    pub const fn from_u64_reduce(v: u64) -> GF257 {
        // R^2 = 2**16 + 2**9 + 1 mod p.
        // With v < 2^64, we have R^2*v < 2^80 + 2^73 + 2^1, which is in
        // range of montyred().
        GF257(GF257::montyred((v as u128) * (GF257::R2 as u128)))
    }

    /// Get the element as an integer, normalized in the 0..p-1
    /// range.
    #[inline(always)]
    pub const fn to_u64(self) -> u64 {
        // Conversion back to normal representation is only a matter of
        // dividing by 2^64 modulo p, and that is exactly what montyred()
        // computes.
        GF257::montyred(self.0 as u128)
    }

    /// Addition in GF(p)
    #[inline(always)]
    const fn add(self, rhs: Self) -> Self {
        // We compute a + b = a - (p - b).
        let (x1, c1) = self.0.overflowing_sub(GF257::MOD - rhs.0);
        let adj = ((c1 as u64) << 8) + (c1 as u64);
        GF257(x1.wrapping_sub(adj as u64))
    }

    /// Subtraction in GF(p)
    #[inline(always)]
    const fn sub(self, rhs: Self) -> Self {
        // See montyred() for details on the subtraction.
        let (x1, c1) = self.0.overflowing_sub(rhs.0);
        let adj = ((c1 as u64) << 8) + (c1 as u64);
        GF257(x1.wrapping_sub(adj))
    }

    /// Negation in GF(p)
    #[inline(always)]
    const fn neg(self) -> Self {
        GF257::ZERO.sub(self)
    }

    /// Halving in GF(p) (division by 2).
    #[inline(always)]
    pub const fn half(self) -> Self {
        // If x is even, then this returned x/2.
        // If x is odd, then this returns (x-1)/2 + (p+1)/2 = (x+p)/2.
        GF257((self.0 >> 1).wrapping_add((self.0 & 1).wrapping_neg() & 0x7FFFFFFFFFFFFF80))
    }

    /// Doubling in GF(p) (multiplication by 2).
    #[inline(always)]
    pub const fn double(self) -> Self {
        self.add(self)
    }

    /// Multiplication in GF(p) by a small integer (less than 2^31).
    #[inline(always)]
    pub const fn mul_small(self, rhs: u32) -> Self {
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
        let (r, c) = xl.overflowing_sub(GF257::MOD - ((xh << 8) + xh));
        let adj = ((c as u64) << 8) + (c as u64);
        GF257(r.wrapping_sub(adj))
    }

    /// Multiplication in GF(p)
    #[inline(always)]
    const fn mul(self, rhs: Self) -> Self {
        // If x < p and y < p, then x*y <= (p-1)^2, and is thus in
        // range of montyred().
        GF257(GF257::montyred((self.0 as u128) * (rhs.0 as u128)))
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
        self.equals(&GF257::ONE)
    }

    /// Test of equality with minus one; return value is 0xFFFFFFFF
    /// if this value is equal to -1 mod p, or 0 otherwise.
    #[inline(always)]
    pub const fn isminusone(self) -> u32 {
        self.equals(&GF257::MINUS_ONE)
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

    /// Legendre symbol: return x^((p-1)/2) (as a GF(p) element).
    pub fn legendre_gfp(self) -> GF257 {
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
        let c1 = l.equals(&GF257::ONE);
        let c2 = l.equals(&GF257::MINUS_ONE);
        let cc1 = (c1 & 1) as i32;
        let cc2 = (c2 & 1) as i32;

        cc1 - cc2
    }

    /// Square root in GF(p); returns (r, cc):
    ///  - If the input is a square, r = sqrt(self) and cc = 0xFFFFFFFF
    ///  - If the input is not a square, r = zero and cc = 0
    /// Which of the two square roots is obtained is unspecified.
    pub fn sqrt(self) -> (Self, u32) {
        // We use that p = 3 mod 4 to compute the square root by
        // x^((p+1)/4) mod p. We have (p+1)/4 = 2^62 - 2^6
        // In the instructions below, we call 'xj' the value x^(2^j-1).

        let x = self;

        // y <- u^((q+1)/4)
        let x2 = x * x.square();
        let x4 = x2 * x2.msquare(2);
        let x6 = x2 * x4.msquare(2);
        let x7 = x * x6.square();
        let x14 = x7 * x7.msquare(7);
        let x28 = x14 * x14.msquare(14);
        let x56 = x28 * x28.msquare(28);
        let y = x56.msquare(6);

        let r = y.square().equals(&x);
        let m = (r as u64) | ((r as u64) << 32);
        (GF257(y.0 & m), m as u32)
    }

    // Conditionally copy the provided value ('a') into self:
    //  - If ctl == 0xFFFFFFFF, then the value of 'a' is copied into self.
    //  - If ctl == 0, then the value of self is unchanged.
    // ctl MUST be equal to 0 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn set_cond(&mut self, a: &Self, ctl: u32) {
        let cw = ((ctl as i32) as i64) as u64;
        self.0 ^= cw & (self.0 ^ a.0);
    }

    /// Negate this value if ctl is 0xFFFFFFFF; leave it unchanged if
    /// ctl is 0x00000000.
    /// The value of ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn set_condneg(&mut self, ctl: u32) {
        let v = -(self as &Self);
        self.set_cond(&v, ctl);
    }

    /// Select a value: this function returns x0 if c == 0, or x1 if
    /// c == 0xFFFFFFFF.
    #[inline(always)]
    pub fn select(x0: GF257, x1: GF257, ctl: u32) -> GF257 {
        let c = (ctl as u64) | ((ctl as u64) << 32);
        GF257(x0.0 ^ (c & (x0.0 ^ x1.0)))
    }

    // Conditionally swap two elements: values a and b are exchanged if
    // ctl == 0xFFFFFFFF, or not exchanged if ctl == 0x00000000. Value
    // ctl MUST be either 0x00000000 or 0xFFFFFFFF.
    #[inline(always)]
    pub fn condswap(a: &mut Self, b: &mut Self, ctl: u32) {
        let cw = ((ctl as i32) as i64) as u64;
        let t = cw & (a.0 ^ b.0);
        a.0 ^= t;
        b.0 ^= t;
    }
}

// We implement all the needed traits to allow use of the arithmetic
// operators on GF(p) values.

impl Add<GF257> for GF257 {
    type Output = GF257;

    #[inline(always)]
    fn add(self, other: GF257) -> GF257 {
        GF257::add(self, other)
    }
}

impl AddAssign<GF257> for GF257 {
    #[inline(always)]
    fn add_assign(&mut self, other: GF257) {
        *self = GF257::add(*self, other);
    }
}

impl Sub<GF257> for GF257 {
    type Output = GF257;

    #[inline(always)]
    fn sub(self, other: GF257) -> GF257 {
        GF257::sub(self, other)
    }
}

impl SubAssign<GF257> for GF257 {
    #[inline(always)]
    fn sub_assign(&mut self, other: GF257) {
        *self = GF257::sub(*self, other);
    }
}

impl Neg for GF257 {
    type Output = GF257;

    #[inline(always)]
    fn neg(self) -> GF257 {
        GF257::neg(self)
    }
}

impl Neg for &GF257 {
    type Output = GF257;

    #[inline(always)]
    fn neg(self) -> GF257 {
        GF257::neg(*self)
    }
}

impl Mul<GF257> for GF257 {
    type Output = GF257;

    #[inline(always)]
    fn mul(self, other: GF257) -> GF257 {
        GF257::mul(self, other)
    }
}

impl MulAssign<GF257> for GF257 {
    #[inline(always)]
    fn mul_assign(&mut self, other: GF257) {
        *self = GF257::mul(*self, other);
    }
}

impl Div<GF257> for GF257 {
    type Output = GF257;

    #[inline(always)]
    fn div(self, other: GF257) -> GF257 {
        GF257::div(self, other)
    }
}

impl DivAssign<GF257> for GF257 {
    #[inline(always)]
    fn div_assign(&mut self, other: GF257) {
        *self = GF257::div(*self, other);
    }
}

// ========================================================================
// Unit tests.

#[cfg(test)]
mod tests {
    use super::GF257;

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

    fn check_GF257_eq(a: GF257, r: u128) {
        assert!(a.to_u64() == (r % (GF257::MOD as u128)) as u64);
    }

    fn test_GF257_ops(a: u64, b: u64) {
        let x = GF257::from_u64_reduce(a);
        let y = GF257::from_u64_reduce(b);
        let wa = a as u128;
        let wb = b as u128;

        check_GF257_eq(x + y, wa + wb);
        check_GF257_eq(x - y, (wa + (GF257::MOD as u128) * 2) - wb);
        check_GF257_eq(-y, (GF257::MOD as u128) * 2 - wb);
        check_GF257_eq(x * y, wa * wb);
        check_GF257_eq(x.square(), wa * wa);
        if a == 0 || a == GF257::MOD {
            check_GF257_eq(x.invert(), 0);
        } else {
            check_GF257_eq(x * x.invert(), 1);
        }
        assert!(x.half().mul_small(2).equals(&x) == 0xFFFFFFFF);
        assert!(x.mul_small(3).equals(&(x + x + x)) == 0xFFFFFFFF);
    }

    #[test]
    fn GF257_ops() {
        for i in 0..10 {
            // let v: u64 = (i as u64) + 0xFFFFFFFEFFFFFFFC;
            let v: u64 = (i as u64) + 0xFFFFFFFFFFFFFEFA;
            if i <= 4 {
                let (x, c) = GF257::from_u64(v);
                assert!(c == 0xFFFFFFFF);
                assert!(x.to_u64() == v);
                let y = GF257::from_u64_reduce(v);
                assert!(y.to_u64() == v);
            } else {
                let v2 = v - GF257::MOD;
                let (x, c) = GF257::from_u64(v);
                assert!(c == 0);
                assert!(x.to_u64() == 0);
                let y = GF257::from_u64_reduce(v);
                assert!(y.to_u64() == v2);
            }
        }

        test_GF257_ops(0, 0);
        test_GF257_ops(0, 1);
        test_GF257_ops(1, 0);
        test_GF257_ops(1, 1);
        test_GF257_ops(0, 0xFFFFFFFF);
        test_GF257_ops(0xFFFFFFFF, 0);
        test_GF257_ops(0xFFFFFFFF, 0xFFFFFFFF);
        test_GF257_ops(0, 0xFFFFFFFF00000000);
        test_GF257_ops(0xFFFFFFFF00000000, 0);
        test_GF257_ops(0xFFFFFFFF00000000, 0xFFFFFFFF00000000);
        let mut prng = PRNG(0);
        for _ in 0..10000 {
            let a = prng.next_u64();
            let b = prng.next_u64();
            test_GF257_ops(a, b);
        }
        assert!(GF257::ZERO.legendre_gfp().iszero() == 0xFFFFFFFF);
        let (s0, c0) = GF257::ZERO.sqrt();
        check_GF257_eq(s0, 0);
        assert!(c0 == 0xFFFFFFFF);
        for _ in 0..1000 {
            let x = GF257::from_u64_reduce((prng.next_u64() >> 1) + 1).square();
            assert!(x.legendre_gfp().equals(&GF257::ONE) == 0xFFFFFFFF);
            let (r1, c1) = x.sqrt();
            assert!(r1.square().equals(&x) == 0xFFFFFFFF);
            assert!(c1 == 0xFFFFFFFF);
            let y = x * GF257::from_u64_reduce(7);
            assert!(y.legendre_gfp().equals(&GF257::MINUS_ONE) == 0xFFFFFFFF);
            let (r2, c2) = y.sqrt();
            assert!(r2.iszero() == 0xFFFFFFFF);
            assert!(c2 == 0);
        }
    }
}
