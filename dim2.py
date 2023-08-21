from collections import namedtuple
from dim1 import CGL, ThetaCGL

ThetaNullPointDim2 = namedtuple("ThetaNullPoint_dim_2", "a b c d")

class ThetaCGLDim2(CGL):
    def __init__(self, domain, **kwds):
        super().__init__(domain, chunk=3, **kwds)

    @staticmethod
    def from_elliptic_curves(E1, E2, sqrt_function=None):
        O1=ThetaCGL(E1)
        O2=ThetaCGL(E2)
        a,b = O1.domain
        c,d = O2.domain
        OO = ThetaNullPointDim2(a*c, a*d, b*c, b*d)
        return ThetaCGLDim2(OO, sqrt_function=sqrt_function)

    @staticmethod
    def hadamard(x,y,z,t):
        return (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)

    def radical_2isogeny(self, bits=[0,0,0]):
        """
        Given a level 2-theta null point, compute a 2-isogeneous theta null
        point
        """

        a, b, c, d = self.domain
        aa = a*a # a*a is faster than a**2 in SageMath
        bb = b*b
        cc = c*c
        dd = d*d
        AA, BB, CC, DD = ThetaCGLDim2.hadamard(a,b,c,d)
        AABB = AA * BB
        AACC = AA * CC
        AADD = AA * DD
        AB = self.sqrt(AABB)
        AC = self.sqrt(AACC)
        AD = self.sqrt(AADD)
        
        if bits[0] == 1:
            AB = - AB
        if bits[1] == 1:
            AC = - AC
        if bits[2] == 1:
            AD = - AD

        anew, bnew, cnew, dnew = ThetaCGLDim2.hadamard(AA, AB, AC, AD)
        O1 = ThetaNullPointDim2(anew, bnew, cnew, dnew)
        return O1

    def advance(self, bits=[0,0,0]):
        O1 = self.radical_2isogeny(bits)
        return ThetaCGLDim2(O1, sqrt_function=self.sqrt_function)

    def to_hash(self):
        a, b, c, d = self.domain
        a_inv = 1/a
        return (b * a_inv, c * a_inv, d * a_inv)
