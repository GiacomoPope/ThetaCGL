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
    F = x.parent()
    x0, x1 = x.list()

    delta = x0**2 + x1**2
    n = fourth_Fp(delta)   

    
    # x0^4 - n*x0^2 + (n^2 - t0) / 8 = 0



if __name__ == "__main__":
    p = 79*2**247 - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])

    for _ in range(100):
        x = randint(0, p)
        a = F([x, 0])
        assert sqrt_Fp2(a**2) in [a, -a]

    for _ in range(100):
        y = randint(0, p)
        a = F([0, y])
        assert sqrt_Fp2(a**2) in [a, -a]

    for _ in range(100):
        a = F.random_element()
        assert sqrt_Fp2(a**2) in [a, -a]    