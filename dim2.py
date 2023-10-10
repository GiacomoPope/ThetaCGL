from collections import namedtuple
from cgl import CGL
from dim1 import ThetaCGL
from utilities import canonical_root

ThetaNullPointDim2 = namedtuple("ThetaNullPoint_dim_2", "a b c d")


class ThetaCGLDim2(CGL):
    def __init__(self, domain, sqrt_function=None, chunk=3, **kwds):
        super().__init__(domain, sqrt_function=sqrt_function, chunk=chunk, **kwds)

    @staticmethod
    def from_elliptic_curves(E1, E2, sqrt_function=None):
        O1 = ThetaCGL(E1)
        O2 = ThetaCGL(E2)
        a, b = O1.domain
        c, d = O2.domain
        OO = ThetaNullPointDim2(a * c, a * d, b * c, b * d)
        return ThetaCGLDim2(OO, sqrt_function=sqrt_function)

    @staticmethod
    def hadamard(x, y, z, t):
        return (x + y + z + t, x - y + z - t, x + y - z - t, x - y - z + t)

    def radical_2isogeny(self, bits=[0, 0, 0]):
        """
        Given a level 2-theta null point, compute a 2-isogeneous theta null
        point
        """
        a, b, c, d = self.domain

        aa = a * a  # a*a is faster than a**2 in SageMath
        bb = b * b
        cc = c * c
        dd = d * d

        AA, BB, CC, DD = ThetaCGLDim2.hadamard(aa, bb, cc, dd)

        AABB = AA * BB
        AACC = AA * CC
        AADD = AA * DD
        AB = self.sqrt(AABB)
        AC = self.sqrt(AACC)
        AD = self.sqrt(AADD)

        if bits[0] == 1:
            AB = -AB
        if bits[1] == 1:
            AC = -AC
        if bits[2] == 1:
            AD = -AD

        a_new, b_new, c_new, d_new = ThetaCGLDim2.hadamard(AA, AB, AC, AD)

        return ThetaNullPointDim2(a_new, b_new, c_new, d_new)

    def advance(self, bits=[0, 0, 0]):
        O1 = self.radical_2isogeny(bits)
        return ThetaCGLDim2(O1, sqrt_function=self.sqrt_function)

    def to_hash(self):
        a, b, c, d = self.domain
        a_inv = 1 / a
        return (b * a_inv, c * a_inv, d * a_inv)


class ThetaCGLDim2Radical4(CGL):
    def __init__(
        self,
        domain,
        sqrt_function=None,
        fourth_root_function=None,
        zeta4=None,
        chunk=6,
        **kwds,
    ):
        super().__init__(domain, chunk=chunk, **kwds)

        if zeta4 is None:
            a = self.domain[0]
            zeta4 = a.parent().gen()
        assert zeta4 * zeta4 == -1
        self.zeta4 = zeta4
        self.sqrt_function = sqrt_function
        self.fourth_root_function = fourth_root_function

    @staticmethod
    def hadamard(x, y, z, t):
        return (x + y + z + t, x - y + z - t, x + y - z - t, x - y - z + t)

    @staticmethod
    def from_elliptic_curves(
        E1, E2, sqrt_function=None, fourth_root_function=None, zeta4=None
    ):
        O1 = ThetaCGL(E1)
        O2 = ThetaCGL(E2)
        a, b = O1.domain
        c, d = O2.domain
        OO = ThetaNullPointDim2(a * c, a * d, b * c, b * d)
        return ThetaCGLDim2Radical4(
            OO,
            sqrt_function=sqrt_function,
            fourth_root_function=fourth_root_function,
            zeta4=zeta4,
            chunk=6,
        )

    def radical_4isogeny(self, bits=[0, 0, 0, 0, 0, 0]):
        """
        Given a level 2-theta null point, compute a 4-isogeneous theta null
        point
        """
        a, b, c, d = self.domain

        aa = a * a  # a*a is faster than a**2 in SageMath
        bb = b * b
        cc = c * c
        dd = d * d

        AA, BB, CC, DD = ThetaCGLDim2.hadamard(aa, bb, cc, dd)

        AABB = AA * BB
        AACC = AA * CC
        BBDD = BB * DD
        CCDD = CC * DD

        ABCD = self.sqrt(AABB * CCDD)
        if bits[0] == 1:
            ABCD = -ABCD

        alpha1_4 = 4 * (2 * ABCD + AABB + CCDD)
        alpha2_4 = 4 * (2 * ABCD + AACC + BBDD)

        alpha1 = self.fourth_root(alpha1_4)
        alpha2 = self.fourth_root(alpha2_4)

        if bits[1] == 1:
            alpha1 = -alpha1
        if bits[2] == 1:
            alpha1 *= self.zeta4

        if bits[3] == 1:
            alpha2 = -alpha2
        if bits[4] == 1:
            alpha2 *= self.zeta4

        alpha3_2 = 8 * (CCDD + ABCD)
        alpha3_2 *= (AACC + ABCD) * CCDD * DD + (BBDD + ABCD) * CCDD * CC
        alpha3 = self.sqrt(alpha3_2)

        projective_factor = CCDD * alpha1 * alpha2

        if bits[5] == 1:
            alpha3 = -alpha3

        a_new, b_new, c_new, d_new = ThetaCGLDim2.hadamard(
            2 * a * projective_factor,
            alpha1 * projective_factor,
            alpha2 * projective_factor,
            alpha3,
        )

        return ThetaNullPointDim2(a_new, b_new, c_new, d_new)

    def advance(self, bits=[0, 0, 0, 0, 0, 0]):
        O1 = self.radical_4isogeny(bits=bits)
        return ThetaCGLDim2Radical4(
            O1,
            fourth_root_function=self.fourth_root_function,
            sqrt_function=self.sqrt_function,
            zeta4=self.zeta4,
        )

    def to_hash(self):
        a, b, c, d = self.domain
        a_inv = 1 / a
        return (b * a_inv, c * a_inv, d * a_inv)
