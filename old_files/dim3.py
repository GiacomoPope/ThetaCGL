from collections import namedtuple
from dim1 import CGL, ThetaCGL

ThetaNullPointDim3 = namedtuple("ThetaNullPoint_dim_3", "a b c d e f g h")


class ThetaCGLDim3(CGL):
    def __init__(self, domain, **kwds):
        super().__init__(domain, chunk=6, **kwds)

    @staticmethod
    def from_elliptic_curves(E1, E2, E3, sqrt_function=None):
        O1 = ThetaCGL(E1)
        O2 = ThetaCGL(E2)
        O3 = ThetaCGL(E3)
        a, b = O1.domain
        c, d = O2.domain
        e, f = O3.domain
        OO = ThetaNullPointDim3(
            a * c * e,
            a * d * e,
            b * c * e,
            b * d * e,
            a * c * f,
            a * d * f,
            b * c * f,
            b * d * f,
        )
        return ThetaCGLDim3(OO, sqrt_function=sqrt_function)

    @staticmethod
    def hadamard(x0, x1, x2, x3, x4, x5, x6, x7):
        y0 = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7
        y1 = x0 - x1 + x2 - x3 + x4 - x5 + x6 - x7
        y2 = x0 + x1 - x2 - x3 + x4 + x5 - x6 - x7
        y3 = x0 - x1 - x2 + x3 + x4 - x5 - x6 + x7
        y4 = x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7
        y5 = x0 - x1 + x2 - x3 - x4 + x5 - x6 + x7
        y6 = x0 + x1 - x2 - x3 - x4 - x5 + x6 + x7
        y7 = x0 - x1 - x2 + x3 - x4 + x5 + x6 - x7
        return (y0, y1, y2, y3, y4, y5, y6, y7)

    @staticmethod
    def eval_F2(x000, x001, x010, x011, x100, x101, x110, x111):
        # should be zero on products
        F = (
            (8192)
            * x111
            * x101
            * x100
            * x000
            * (-x000 * x010 - x001 * x011 + x100 * x110 + x101 * x111)
            * (-x000 * x010 + x001 * x011 - x100 * x110 + x101 * x111)
            * (x000 * x010 - x001 * x011 - x100 * x110 + x101 * x111)
            * (x000 * x010 + x001 * x011 + x100 * x110 + x101 * x111)
            * (-x001 * x010 * x101 * x110 + x000 * x011 * x100 * x111)
            * (
                -(x000**3) * x011**3
                + 2 * x000 * x001**2 * x011 * x100**2
                - 2 * x001 * x010**3 * x101 * x111
                + x000 * x011 * x100**2 * x111**2
            )
        )
        return F

    @staticmethod
    def eval_F(x0, x1, x2, x3, x4, x5, x6, x7):
        c0 = x0**2
        c1 = x1**2
        c2 = x2**2
        c3 = x3**2
        c4 = x4**2
        c5 = x5**2
        c6 = x6**2
        c7 = x7**2
        a0, a1, a2, a3, a4, a5, a6, a7 = ThetaCGLDim3.hadamard(
            c0, c1, c2, c3, c4, c5, c6, c7
        )

        b0 = 2 * (x0 * x4 + x1 * x5 + x2 * x6 + x3 * x7)
        b1 = 2 * (x0 * x4 - x1 * x5 + x2 * x6 - x3 * x7)
        b2 = 2 * (x0 * x4 + x1 * x5 - x2 * x6 - x3 * x7)
        b3 = 2 * (x0 * x4 - x1 * x5 - x2 * x6 + x3 * x7)

        R1 = a0 * a1 * a2 * a3
        R2 = b0 * b1 * b2 * b3
        R3 = a4 * a5 * a6 * a7

        F = R1**2 + R2**2 + R3**2 - 2 * (R1 * R2 + R1 * R3 + R2 * R3)

        return F

    @staticmethod
    def last_sqrt(self, c0, c1, c2, c3, c4, c5, c6, c7, x0, x1, x2, x3, x4, x5, x6):
        a0, a1, a2, a3, a4, a5, a6, a7 = ThetaCGLDim3.hadamard(
            c0, c1, c2, c3, c4, c5, c6, c7
        )
        R1 = a0 * a1 * a2 * a3
        R3 = a4 * a5 * a6 * a7

        c04 = c0 * c4
        c15 = c1 * c5
        c26 = c2 * c6
        c37 = c3 * c7
        c0246 = c04 * c26
        c1357 = c15 * c37

        term = 16 * (c04 - c15 + c26 - c37) ** 2 - 64 * (c0246 + c1357)

        num = (R1 + R3 - term) ** 2 + 16384 * c0246 * c1357 - 4 * R1 * R3
        den = 256 * (R1 + R3 - term) * x0 * x1 * x2 * x3 * x4 * x5 * x6

        if den == 0:
            print("method doesn't work, both square-roots are correct")
            x7 = self.sqrt(c7)
        else:
            print("method works")
            x7 = num / den
            assert x7 * x7 == c7, "square-root not computed correctly"

        return x7

    def last_sqrt_slow(
        self, y0, y1, y2, y3, y4, y5, y6, y7, x0, x1, x2, x3, x4, x5, x6
    ):
        # input are the squares of the dual theta nullpoint coordinates
        # and the first seven coordinates of the dual theta nullpoint

        assert y0 == x0**2
        assert y1 == x1**2
        assert y2 == x2**2
        assert y3 == x3**2
        assert y4 == x4**2
        assert y5 == x5**2
        assert y6 == x6**2

        x7 = self.sqrt(y7)

        if not ThetaCGLDim3.eval_F(x0, x1, x2, x3, x4, x5, x6, x7) == 0:
            x7 = -x7
            assert ThetaCGLDim3.eval_F(x0, x1, x2, x3, x4, x5, x6, x7) == 0
        # print("on the boundary?", ThetaCGLDim3.eval_F2(x0,x1,x2,x3,x4,x5,x6,x7))
        if ThetaCGLDim3.eval_F(x0, x1, x2, x3, x4, x5, x6, -x7) == 0:
            print("square-root not uniquely determined")
        else:
            print("unique sqrt")

        assert x7 * x7 == y7, "square-root incorrect"

        return x7

    def radical_2isogeny(self, bits=[0, 0, 0, 0, 0, 0]):
        """
        Given a level 2-theta null point, compute a 2-isogeneous theta null
        point
        """
        a, b, c, d, e, f, g, h = self.domain

        aa = a * a  # a*a is faster than a**2 in SageMath
        bb = b * b
        cc = c * c
        dd = d * d
        ee = e * e
        ff = f * f
        gg = g * g
        hh = h * h
        AA, BB, CC, DD, EE, FF, GG, HH = ThetaCGLDim3.hadamard(
            aa, bb, cc, dd, ee, ff, gg, hh
        )

        AABB = AA * BB
        AACC = AA * CC
        AADD = AA * DD
        AAEE = AA * EE
        AAFF = AA * FF
        AAGG = AA * GG
        AAHH = AA * HH
        AB = self.sqrt(AABB)
        AC = self.sqrt(AACC)
        AD = self.sqrt(AADD)
        AE = self.sqrt(AAEE)
        AF = self.sqrt(AAFF)
        AG = self.sqrt(AAGG)

        if bits[0] == 1:
            AB = -AB
        if bits[1] == 1:
            AC = -AC
        if bits[2] == 1:
            AD = -AD
        if bits[3] == 1:
            AE = -AE
        if bits[4] == 1:
            AF = -AF
        if bits[5] == 1:
            AG = -AG

        AH = ThetaCGLDim3.last_sqrt(
            self,
            AA * AA,
            AABB,
            AACC,
            AADD,
            AAEE,
            AAFF,
            AAGG,
            AAHH,
            AA,
            AB,
            AC,
            AD,
            AE,
            AF,
            AG,
        )

        anew, bnew, cnew, dnew, enew, fnew, gnew, hnew = ThetaCGLDim3.hadamard(
            AA, AB, AC, AD, AE, AF, AG, AH
        )
        O1 = ThetaNullPointDim3(anew, bnew, cnew, dnew, enew, fnew, gnew, hnew)
        return O1

    def advance(self, bits=[0, 0, 0, 0, 0, 0]):
        O1 = self.radical_2isogeny(bits)
        return ThetaCGLDim3(O1, sqrt_function=self.sqrt_function)

    def to_hash(self):
        a, b, c, d, e, f, g, h = self.domain
        a_inv = 1 / a
        return (
            b * a_inv,
            c * a_inv,
            d * a_inv,
            e * a_inv,
            f * a_inv,
            g * a_inv,
            h * a_inv,
        )
