from sage.all import GF, proof

proof.all(False)


def invert_or_zero(x):
    if x == 0:
        return 0
    return 1 / x


def sqrt_Fp(x):
    p = x.parent().characteristic()
    exp = (p + 1) // 4

    r = x**exp
    if r * r != x:
        return 0
    if int(r) % 2 != 0:
        return -r
    return r


def sqrt_Fp2(x):
    F = x.parent()
    x0, x1 = x.list()

    if x1 == 0:
        lx0 = x0.is_square()
        if lx0:
            y0 = sqrt_Fp(x0)
            return F([y0, 0])
        else:
            y1 = sqrt_Fp(-x0)
            return F([0, y1])

    delta = x0**2 + x1**2
    sqrt_delta = sqrt_Fp(delta)

    y02 = (x0 + sqrt_delta) / 2
    if not y02.is_square():
        y02 -= sqrt_delta

    y0 = sqrt_Fp(y02)
    y1 = x1 / (y0 + y0)

    return F([y0, y1])


def fourth_Fp(x):
    p = x.parent().characteristic()
    exp = (p + 1) // 8

    r = x**exp
    if r * r * r * r != x:
        return 0

    if int(r) % 2 != 0:
        return -r
    return r


def fourth_Fp2(x):
    """
    The goal is to find elements y0, y1 in Fp
    such that x = x0 + ix1 = (y0 + iy1)^4

    As the norm is multiplicative, we have that
    (x0^2 + x1^2) = (y0^2 + y1^2)^4 and so with
    one fourth root in Fp we have (y0^2 + y1^2)

    Expanding out the fourth power, we have

    x0 = y0^4 - 6y0^2y1^2 + y1^4
    x1 = 4y0y1(y0^2 - y1^2)

    Which together with

    n = y0^2 + y1^2

    Allows us to find a quadratic equation in y0^2

    8y0^4 - 8ny0^2 + n^2 - x0 = 0

    We can find the roots of this with on square-root
    in Fp to find

    y0^2 = [8n + sqrt(32(n^2 + x0))] / 16

    y0 is recovered by one last sqrt in Fp and y1 from
    an inversion

    y1 = x1 / 4y0(2y0^2 - n)
    """
    F = x.parent()
    x0, x1 = x.list()

    delta = x0**2 + x1**2
    n = fourth_Fp(delta)
    assert n * n * n * n == delta, "Fourth root didnt work"

    disc = (n**2 + x0) / 2
    disc_sqrt = sqrt_Fp(disc)
    assert disc_sqrt * disc_sqrt == disc, "disc_sqrt didnt work"

    # We do not know which of n or -n
    # is correct, test with legendre
    y02 = (n + disc_sqrt) / 2

    if not y02.is_square():
        y02 -= n
        n = -n

    if y02 == 0:
        y0 = sqrt_Fp(n)
        assert y0 * y0 == n, "n sqrt didnt work"
    else:
        y0 = sqrt_Fp(y02)
        assert y0 * y0 == y02, "y0^2 sqrt didnt work"

    y1 = x1 * invert_or_zero(4 * y0 * disc_sqrt)

    # Handle case with x1 = 0
    if x1 == 0 and disc == 0:
        y1 = y0

    r = F([y0, y1])
    assert r**4 == x, "Whole thing is broken??"
    return r

def eighth_Fp(x):
    p = x.parent().characteristic()
    exp = (p + 1) // 16

    r = x**exp
    if r**8 != x:
        raise ValueError("eighth root failed in Fp")

    if int(r) % 2 != 0:
        return -r
    return r

def eighth_Fp2(x):
    """
    Compute 8th roots in Fp2
    """
    def inner_sqrt(A0, B0, N1, check_square=False):
        
        C1 = (A0 + N1) / 2

        # For the last sqrt we need to ensure C1 is a square
        if check_square and not C1.is_square():
            C1 = (A0 - N1) / 2
            N1 = -N1
        assert C1.is_square(), "C1 must be a square..."
        
        # When C1 is zero only when x1 is zero
        if C1.is_zero():
            assert x1 == 0
            A1 = sqrt_Fp(N1)
            assert A1 * A1 == N1, "N1 sqrt failed..."
        else:
            A1 = sqrt_Fp(C1)
            assert A1 * A1 == C1, "A1 sqrt failed..."

        B1 = B0 * invert_or_zero(2*A1)
        
        # In this case, we must switch A1, B1
        if C1.is_zero():
            assert x1 == 0
            A1, B1 = B1, A1

        return A1, B1

    F = x.parent()
    x0, x1 = x.list()
    delta = x0**2 + x1**2
    n = eighth_Fp(delta)
    assert n**8 == delta, "Eighth root didnt work"

    # initialization
    A0 = x0
    B0 = x1

    A1, B1 = inner_sqrt(A0, B0, n**4)    
    A2, B2 = inner_sqrt(A1, B1, n**2)
    A3, B3 = inner_sqrt(A2, B2, n**1, check_square=True)

    return F([A3, B3])


if __name__ == "__main__":
    p = 127
    assert p % 16 == 15
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    i = F.gen()

    for b0 in range(p):
        for b1 in range(p):
            b = F([b0, b1])
            a = b**8
            try:
                x = eighth_Fp2(a)
            except Exception as e:
                print(f"{e = }")
                print(f"{a = }")
                print(f"{b = }")
                continue
            assert x**8 == a, f"{x = }, {a = }, {a.nth_root(8, all=True) = }"

    # Big prime
    p = 5*2**248 - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    for _ in range(1000):
        b = F.random_element()
        a = b**8
        x = eighth_Fp2(a)
        assert x**8 == a
