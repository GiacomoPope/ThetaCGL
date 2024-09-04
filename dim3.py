from collections import namedtuple
from dim1 import CGL, ThetaCGL
from sage.all import vector

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
        # We changed the ordering to be compatible with our notes.
        # OO = ThetaNullPointDim3(a*c*e, a*d*e, b*c*e, b*d*e, a*c*f, a*d*f, b*c*f, b*d*f)
        OO = ThetaNullPointDim3(
            a * c * e,
            b * c * e,
            a * d * e,
            b * d * e,
            a * c * f,
            b * c * f,
            a * d * f,
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
    def last_sqrt(
        self, x0, x1, x2, x3, x4, x5, x6, x7, y0, y1, y2, y3, y4, y5, y6
    ):
        b0, b1, b2, b3, b4, b5, b6, b7 = ThetaCGLDim3.hadamard(
            x0, x1, x2, x3, x4, x5, x6, x7
        )
        R1 = b0 * b1 * b2 * b3
        R3 = b4 * b5 * b6 * b7

        x04 = x0 * x4
        x15 = x1 * x5
        x26 = x2 * x6
        x37 = x3 * x7
        x0246 = x04 * x26
        x1357 = x15 * x37

        t0 = 16 * (x04 - x15 + x26 - x37) ** 2 - 64 * (x0246 + x1357)
        T = R1 + R3 - t0
        y = y0 * y1 * y2 * y3 * y4 * y5 * y6

        if T != 0:
            t1 = T**2 + 16384 * x0246 * x1357 - 4 * R1 * R3
            t2 = 256 * T * y
        else:
            a0, a1, a2, a3, a4, a5, a6, a7 = self.domain
            t1 = -64 * a0 * a1 * a2 * a3 * a4 * a5 * a6 * a7 * x0**2
            t2 = y

        y0, y1, y2, y3, y4, y5, y6 = [t2 * y for y in [y0, y1, y2, y3, y4, y5, y6]]
        y7 = t1

        return y0, y1, y2, y3, y4, y5, y6, y7

    def squared_thetas(self):
        converter = {
            0: vector([0, 0, 0]),
            1: vector([1, 0, 0]),
            2: vector([0, 1, 0]),
            3: vector([1, 1, 0]),
            4: vector([0, 0, 1]),
            5: vector([1, 0, 1]),
            6: vector([0, 1, 1]),
            7: vector([1, 1, 1]),
        }

        a, b, c, d, e, f, g, h = self.domain
        Theta_list = [a, b, c, d, e, f, g, h]

        squared_thetas_list = []

        for chi in range(0, 8):
            for i in range(0, 8):
                if (converter[chi] * converter[i]) % 2 == 0:
                    U_chi_i = 0

                    for t in range(0, 8):
                        chi_t = (-1) ** (converter[chi] * converter[t])
                        vec = (converter[i] + converter[t]) % 2
                        index = vec[0] + 2 * vec[1] + 4 * vec[2]
                        U_chi_i += chi_t * Theta_list[index] * Theta_list[t]
                    # if U_chi_i == 0:
                    #    print(chi, i)
                    squared_thetas_list.append(U_chi_i)
                    # print(80*"=")

        return squared_thetas_list

    def label(self):
        # Type1: plane quartic, Type2: hyperelliptic, Type3:product of dimension 1 and 2, Type4: product of elliptic curves
        squared_thetas_list = self.squared_thetas()
        zeros = 0
        for el in squared_thetas_list:
            if el == 0:
                zeros += 1
        if zeros == 0:
            return "Plane quartic"
        elif zeros == 1:
            return "Hyperelliptic curve"
        elif zeros == 6:
            return "Product of a hyperelliptic curve and an elliptic curve"
        elif zeros == 9:
            return "Product of three elliptic curves"
        else:
            print("unexpected case: number of zeros is", zeros)
            return None

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

        lam = 0
        for X in [AA, BB, CC, DD, EE, FF, GG, HH]:
            if X != 0:
                lam = X
                break
        if lam == 0:
            raise ValueError("All coordinates are zero")

        AAAA = lam * AA
        AABB = lam * BB
        AACC = lam * CC
        AADD = lam * DD
        AAEE = lam * EE
        AAFF = lam * FF
        AAGG = lam * GG
        AAHH = lam * HH

        AB = self.sqrt(AABB)
        AC = self.sqrt(AACC)
        AD = self.sqrt(AADD)
        AE = self.sqrt(AAEE)
        AF = self.sqrt(AAFF)
        AG = self.sqrt(AAGG)

        xx = [AB, AC, AD, AE, AF, AG]
        zero_indices = [i for i in range(6) if xx[i] == 0]

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

        # print(f"{AA == lam = }")

        if zero_indices:
            # in this case, we can still choose the last square-root (?)
            print("we have redundant bits at ", zero_indices)
            AH = self.sqrt(AAHH)
            if bits[zero_indices[0]] == 1:
                AH = -AH
        else:
            AA, AB, AC, AD, AE, AF, AG, AH = ThetaCGLDim3.last_sqrt(
                self,
                AAAA,
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

    def eval_F2(self, x000, x001, x010, x011, x100, x101, x110, x111):
        # should be zero on products
        f1 = x111 * x101 * x100 * x000
        f2 = -x000 * x010 - x001 * x011 + x100 * x110 + x101 * x111
        f3 = -x000 * x010 + x001 * x011 - x100 * x110 + x101 * x111
        f4 = x000 * x010 - x001 * x011 - x100 * x110 + x101 * x111
        f5 = x000 * x010 + x001 * x011 + x100 * x110 + x101 * x111
        f6 = -x001 * x010 * x101 * x110 + x000 * x011 * x100 * x111
        f7 = (
            -(x000**3) * x011**3
            + 2 * x000 * x001**2 * x011 * x100**2
            - 2 * x001 * x010**3 * x101 * x111
            + x000 * x011 * x100**2 * x111**2
        )
        return [x000, x001, x010, x011, x100, x101, x110, x111, f2, f3, f4, f5, f6, f7]

    def advance(self, bits=[0, 0, 0, 0, 0, 0]):
        O1 = self.radical_2isogeny(bits)
        O1 = ThetaCGLDim3(O1, sqrt_function=self.sqrt_function)
        return O1

    def bit_string(self, message):
        r = self
        for x in range(0, len(message), self.chunk):
            bits = message[x : x + self.chunk]
            if len(bits) < self.chunk:
                # pad bits if too short
                bits += [0] * (self.chunk - len(bits))
            try:
                r = r.advance(bits)
            except Exception as e:
                raise ValueError(f"Something went wrong with: {bits}")
        return r

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
