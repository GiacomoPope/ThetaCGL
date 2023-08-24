from sage.all import GF, proof, randint

proof.all(False)

def sqrt_Fp(x):
    p = x.parent().characteristic()
    exp = (p+1) // 4

    assert x.is_square()
    return x ** exp

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
    exp = (p+1) // 8

    return x ** exp

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

    n = y0^2 - y1^2

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

    disc = (n**2 + x0) / 2
    disc_sqrt = sqrt_Fp(disc)

    # We do not know which of n or -n
    # is correct, test with legendre
    y02 = ( n + disc_sqrt) / 2

    if not y02.is_square():
        y02 -= n
        n = -n

    y0 = sqrt_Fp(y02)

    # When we have
    # (y02 + y02 - n) = 0
    # Then we have y0^2 = y1^2
    gamma = y02 + y02 - n
    if gamma.is_zero():
        return F([y0, y0])

    y1 = x1 / (4*y0*gamma)
    return F([y0, y1])

if __name__ == "__main__":
    # p = 79*2**247 - 1

    p = 167
    assert p % 8 == 7
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    i = F.gen()

    for b0 in range(p):
        for b1 in range(p):
            b = F([b0, b1])
            a = b**4

            try:
                x = fourth_Fp2(a)
            except Exception as e:
                print(e)
                print(f"{a = }")
                pass
            assert x**4 == a, f"{a = }, {b = }"