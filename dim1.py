from collections import namedtuple
from sage.all import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

ThetaNullPoint = namedtuple("ThetaNullPoint_dim_1", "a b")

class CGL:
    def __init__(self, domain, sqrt_function=None):
        self.domain=domain
        self.sqrt_function=sqrt_function

    def __repr__(self):
        return f"CGL with domain={self.domain}"

    def sqrt(self, x):
        if self.sqrt_function is None:
            return(x.sqrt())
        else:
            return(self.sqrt_function(x))

    def advance(self, bit=0):
        pass

    def bit_string(self, bits):
        r=self
        for bit in bits:
            r=r.advance(bit)
        return r

    def to_hash():
        pass

    def hash(self, bits):
        r=self.bit_string(bits)
        return r.to_hash()

def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError("The elliptic curve E is not in the Montgomery model.")
    return A


class ThetaCGL(CGL):
    def __init__(self, domain, sqrt_function=None):
        super().__init__(domain, sqrt_function=sqrt_function)
        if isinstance(self.domain, EllipticCurve_generic):
            self.domain=self.montgomery_curve_to_theta_null_point(domain)

    def montgomery_curve_to_theta_null_point(self, E):
        """
        From an elliptic curve in Montgomery form, compute a theta null point
        (a:b).
        Let T1=(1:1) the canonical point of 4-torsion in Montgomery coordinates
        and T2 such that (T1,T2) forms a symplectic basis of E[4].
        There are 4 choices of T2, giving 4 different theta null points.
        """
        # Extract A from curve equation
        A = montgomery_coefficient(E)
    
        # alpha is a root of
        # x^2 + Ax + 1
        disc = A * A - 4
        assert disc.is_square()
    
        d = self.sqrt(A * A - 4)
        alpha = (-A + d) / 2
    
        # The theta coordinates a^2 and b^2
        # are related to alpha
        aa = alpha + 1
        bb = alpha - 1
    
        aabb = aa * bb
        ab = self.sqrt(aabb)
    
        # We aren't given (a,b) rational, but
        # (a/b) is rational so we use
        # (ab : b^2) as the theta null point
        O0 = ThetaNullPoint(ab, bb)
        return O0

    def radical_2isogeny(self, sign=1):
        """
        Given a level 2-theta null point, compute a 2-isogeneous theta null
        point
        """

        a,b=self.domain
        aa=a**2
        bb=b**2
        AA=aa+bb
        BB=aa-bb
        AABB=AA*BB
        AB=self.sqrt(AABB)
        AB=sign*AB
        anew=AA+AB
        bnew=AA-AB
        O1=ThetaNullPoint(anew, bnew)
        return O1

    def advance(self, bit=0):
        sign=1 if bit == 1 else -1
        O1=self.radical_2isogeny(sign)
        return ThetaCGL(O1, sqrt_function=self.sqrt_function)

    def to_hash(self):
        a,b=self.domain
        return b/a
