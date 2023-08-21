from collections import namedtuple

from sage.all import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

from utilities import montgomery_coefficient

ThetaNullPoint = namedtuple("ThetaNullPoint_dim_1", "a b")


class CGL:
    def __init__(self, domain, sqrt_function=None, chunk=1):
        self.domain = domain
        self.sqrt_function = sqrt_function
        self.chunk = chunk

    def __repr__(self):
        return f"CGL with domain={self.domain}"

    def sqrt(self, x):
        if self.sqrt_function is None:
            r = x.sqrt()
        else:
            r = self.sqrt_function(x)

        # TODO:
        # This is stupid, but making the sqrt canonical 
        # helps with the comparison
        return min(r, -r)

    def advance(self, bits=None):
        pass

    def bit_string(self, message):
        r = self
        for x in range(0, len(message), self.chunk):
            bits=message[x:x+self.chunk]
            if len(bits) < self.chunk:
                # pad bits if too short
                bits += [0]*(self.chunk - len(bits))
            r = r.advance(bits)
        return r

    def to_hash():
        pass

    def hash(self, bits):
        r = self.bit_string(bits)
        return r.to_hash()


class ThetaCGL(CGL):
    def __init__(self, domain, **kwds):
        super().__init__(domain, **kwds)
        if isinstance(self.domain, EllipticCurve_generic):
            self.domain = self.montgomery_curve_to_theta_null_point(domain)

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

        d = self.sqrt(disc)
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

    def to_montgomery_curve(self):
        """
        Given a level 2 theta null point (a:b), compute a Montgomery curve equation.
        We use the model where the 4-torsion point (1:0) above (a:-b) is sent
        to (1:1) in Montgomery coordinates.
        """

        a,b = self.domain
        
        aa = a**2
        bb = b**2

        T1 = aa + bb
        T2 = aa - bb

        # Montgomery coefficient
        A = -(T1**2 + T2**2) / (T1 * T2)

        # Construct curve
        F = b.parent()
        E = EllipticCurve(F, [0, A, 0, 1, 0])
        return E

    def j_invariant(self):
        """
        Give the j-invariant of our current theta null point
        """

        return self.to_montgomery_curve().j_invariant()

    @staticmethod
    def hadamard(x,z):
        return (x+z, x-z)

    def radical_2isogeny(self, bits=[0]):
        """
        Given a level 2-theta null point, compute a 2-isogeneous theta null
        point
        """

        a, b = self.domain
        aa = a*a # a*a is faster than a**2 in SageMath
        bb = b*b
        AA,BB = ThetaCGL.hadamard(aa, bb)
        AABB = AA * BB
        AB = self.sqrt(AABB)
        if bits[0] == 1:
            AB = - AB
        anew, bnew = ThetaCGL.hadamard(AA, AB)
        O1 = ThetaNullPoint(anew, bnew)
        return O1

    def advance(self, bits=[0]):
        O1 = self.radical_2isogeny(bits=bits)
        return ThetaCGL(O1, sqrt_function=self.sqrt_function)

    def to_hash(self):
        a, b = self.domain
        return b / a

# Experimental formula for 4-radical isogeny
class ThetaCGLRadical4(ThetaCGL):
    def __init__(self, domain, sqrt4_function = None, zeta = None, **kwds):
        super().__init__(domain, chunk=2, **kwds)
        if isinstance(self.domain, EllipticCurve_generic):
            self.domain = self.montgomery_curve_to_theta_null_point(domain)

        if zeta is None:
            a,b=self.domain
            zeta = a.base_ring().gen()
        assert zeta * zeta == -1
        self.zeta=zeta
        self.sqrt4_function = sqrt4_function

    # Here sqrt is the fourth root power
    def sqrt4(self, x):
        if self.sqrt4_function is None:
            r = self.sqrt(self.sqrt(x))
        else:
            r = self.sqrt4_function(x)

        zeta = self.zeta
        # normalize the fourth root
        return min(r, zeta*r, -r, -zeta*r)

    def radical_4isogeny(self, bits=[0,0]):
        """
        Given a level 2-theta null point, compute a 4-isogeneous theta null
        point
        """

        a, b = self.domain
        aa = a*a # a*a is faster than a**2 in SageMath
        bb = b*b
        AA,BB = ThetaCGL.hadamard(aa, bb)
        AABB = AA * BB
        factor = self.sqrt4(AABB) #fourth root
        if bits[0] == 1:
            factor = - factor
        if bits[1] == 1:
            factor = self.zeta * factor
        anew=a+b+factor
        bnew=a+b-factor
        anew, bnew = ThetaCGL.hadamard(anew, bnew) # I think we need an hadamard?
        O1 = ThetaNullPoint(anew, bnew)
        return O1

    def advance(self, bits=[0, 0]):
        O1 = self.radical_4isogeny(bits=bits)
        return ThetaCGL(O1, sqrt_function=self.sqrt_function)
