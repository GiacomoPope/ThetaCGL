use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use core::convert::TryFrom;

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

    /// GF(p) modulus: p = 2^64 - 2^32 + 1
    /// GF(p) modulus: p = 2^64 - 257
    // pub const MOD: u64 = 0xFFFFFFFF00000001;
    // TODO: take C into account
    pub const MOD: u64 = 0xfffffffffffffeff;

    /// Element 0 in GF(p).
    pub const ZERO: GFp = GFp::from_u64_reduce(0);

    /// Element 1 in GF(p).
    pub const ONE: GFp = GFp::from_u64_reduce(1);

    /// Element -1 in GF(p).
    pub const MINUS_ONE: GFp = GFp::from_u64_reduce(GFp::MOD - 1);

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
    pub fn legendre(self) -> GFp {
        // (p-1)/2 = 0x7FFFFFFF80000000
        // TODO:
        // (p-1)/2 = 0x7fffffffffffff7f
        let x = self;
        let x2 = x * x.square();
        let x4 = x2 * x2.msquare(2);
        let x8 = x4 * x4.msquare(4);
        let x16 = x8 * x8.msquare(8);
        let x32 = x16 * x16.msquare(16);
        x32.msquare(31)
    }

    // For g = 7^(2^32-1) mod p = 1753635133440165772 (a primitive 2^32 root
    // of 1 in GF(p)), we precompute GG[i] = g^(2^i) for i = 0 to 31.
    const GG: [GFp; 32] = [
        GFp::from_u64_reduce( 1753635133440165772),
        GFp::from_u64_reduce( 4614640910117430873),
        GFp::from_u64_reduce( 9123114210336311365),
        GFp::from_u64_reduce(16116352524544190054),
        GFp::from_u64_reduce( 6414415596519834757),
        GFp::from_u64_reduce( 1213594585890690845),
        GFp::from_u64_reduce(17096174751763063430),
        GFp::from_u64_reduce( 5456943929260765144),
        GFp::from_u64_reduce( 9713644485405565297),
        GFp::from_u64_reduce(16905767614792059275),
        GFp::from_u64_reduce( 5416168637041100469),
        GFp::from_u64_reduce(17654865857378133588),
        GFp::from_u64_reduce( 3511170319078647661),
        GFp::from_u64_reduce(18146160046829613826),
        GFp::from_u64_reduce( 9306717745644682924),
        GFp::from_u64_reduce(12380578893860276750),
        GFp::from_u64_reduce( 6115771955107415310),
        GFp::from_u64_reduce(17776499369601055404),
        GFp::from_u64_reduce(16207902636198568418),
        GFp::from_u64_reduce( 1532612707718625687),
        GFp::from_u64_reduce(17492915097719143606),
        GFp::from_u64_reduce(  455906449640507599),
        GFp::from_u64_reduce(11353340290879379826),
        GFp::from_u64_reduce( 1803076106186727246),
        GFp::from_u64_reduce(13797081185216407910),
        GFp::from_u64_reduce(17870292113338400769),
        GFp::from_u64_reduce(        549755813888),
        GFp::from_u64_reduce(      70368744161280),
        GFp::from_u64_reduce(17293822564807737345),
        GFp::from_u64_reduce(18446744069397807105),
        GFp::from_u64_reduce(     281474976710656),
        GFp::from_u64_reduce(18446744069414584320)
    ];

    /// Square root in GF(p); returns (r, cc):
    ///  - If the input is a square, r = sqrt(self) and cc = 0xFFFFFFFFFFFFFFFF
    ///  - If the input is not a square, r = zero and cc = 0
    /// Which of the two square roots is obtained is unspecified.
    pub fn sqrt(self) -> (Self, u64) {
        // We use a constant-time Tonelli-Shanks. 
        // Input: x
        // Output: (sqrt(x), -1) if x is QR, or (0, 0) otherwise
        // Definitions:
        //    modulus: p = q*2^n + 1 with q odd (here, q = 2^32 - 1 and n = 32)
        //    g is a primitive 2^n root of 1 in GF(p)
        //    GG[j] = g^(2^j)  (j = 0 to n-1, precomputed)
        // Init:
        //    r <- x^((q+1)/2)
        //    v <- x^q
        // Process:
        //    for i = n-1 down to 1:
        //        w = v^(2^(i-1))   (with i-1 squarings)
        //        if w == -1 then:
        //            v <- v*GG[n-i]
        //            r <- r*GG[n-i-1]
        //    if v == 0 or 1 then:
        //        return (r, -1)
        //    else:
        //        return (0, 0)  (no square root)

        let x = self;

        // r <- u^((q+1)/2)
        // v <- u^q
        let x2 = x * x.square();
        let x4 = x2 * x2.msquare(2);
        let x5 = x * x4.square();
        let x10 = x5 * x5.msquare(5);
        let x15 = x5 * x10.msquare(5);
        let x16 = x * x15.square();
        let x31 = x15 * x16.msquare(15);
        let mut r = x * x31;
        let mut v = x * x31.square();

        for i in (1..32).rev() {
            let w = v.msquare((i - 1) as u32);
            let cc = w.equals(GFp::MINUS_ONE);
            v = GFp(v.0 ^ (cc & (v.0 ^ (v * GFp::GG[32 - i]).0)));
            r = GFp(r.0 ^ (cc & (r.0 ^ (r * GFp::GG[31 - i]).0)));
        }
        let m = v.iszero() | v.equals(GFp::ONE);
        (GFp(r.0 & m), m)
    }

    /// Select a value: this function returns x0 if c == 0, or x1 if
    /// c == 0xFFFFFFFFFFFFFFFF.
    #[inline(always)]
    pub fn select(c: u64, x0: GFp, x1: GFp) -> GFp {
        GFp(x0.0 ^ (c & (x0.0 ^ x1.0)))
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
            println!("{}", i);
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
        assert!(GFp::ZERO.legendre().iszero() == 0xFFFFFFFFFFFFFFFF);
        let (s0, c0) = GFp::ZERO.sqrt();
        check_gfp_eq(s0, 0);
        assert!(c0 == 0xFFFFFFFFFFFFFFFF);
        for _ in 0..1000 {
            let x = GFp::from_u64_reduce((prng.next_u64() >> 1) + 1).square();
            assert!(x.legendre().equals(GFp::ONE) == 0xFFFFFFFFFFFFFFFF);
            let (r1, c1) = x.sqrt();
            assert!(r1.square().equals(x) == 0xFFFFFFFFFFFFFFFF);
            assert!(c1 == 0xFFFFFFFFFFFFFFFF);
            let y = x * GFp::from_u64_reduce(7);
            assert!(y.legendre().equals(GFp::MINUS_ONE) == 0xFFFFFFFFFFFFFFFF);
            let (r2, c2) = y.sqrt();
            assert!(r2.iszero() == 0xFFFFFFFFFFFFFFFF);
            assert!(c2 == 0);
        }
    }

}