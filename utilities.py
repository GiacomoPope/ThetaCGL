# ========================== #
#     Montgomery Helpers     #
# ========================== #

def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError("The elliptic curve E is not in the Montgomery model.")
    return A

# ============================================ #
#     Fast square root and quadratic roots     #
# ============================================ #

def canonical_root(a):
    """
    Very stupid and slow way, but it
    makes the sqrt match rust for all
    cases
    """
    a0, a1 = a.list()
    if a0.is_zero() and (int(a1) % 2) == 1:
        return -a
    if (int(a0) % 2) == 1:
        return -a
    return a

def new_sqrt_Fp(a, exp=None):
    """
    Shank's algorithm for sqrt
    Alg 2 in ia.cr/2012/685
    """
    if exp is None:
        p = a.base().characteristic()
        exp = (p - 3) // 4

    a1 = a**exp
    return a * a1

def new_sqrt_Fp_irrational(a, exp=None):
    """
    When we have an element a in Fp2 which 
    is of the form a + i*0 then the root of
    a in Fp may be irrational, so we need
    to account for this

    Modification of alg 2 in ia.cr/2012/685
    """
    if exp is None:
        p = a.base().characteristic()
        exp = (p - 3) // 4

    a1 = a**exp
    x = a * a1
    a0 = a1 * x
    if a0 == -1:
        return 0, x
    return x, 0

def new_sqrt_Fp2(a):
    """
    Complex method for roots

    Algorithm 8 in ia.cr/2012/685
    """
    F = a.parent()
    p = F.characteristic()
    exp = (p - 3) // 4

    assert a.is_square(), "Bad from the beginning!"

    # Extract out Fp coeffs
    a0, a1 = a.list()

    # Easy case, we're already in Fp
    if a1 == 0:
        x0, x1 = new_sqrt_Fp_irrational(a0, exp=exp)
        return F([x0, x1])

    # Otherwise
    a0a0 = a0*a0
    a1a1 = a1*a1

    alpha = a0a0 + a1a1

    alpha = new_sqrt_Fp(alpha, exp=exp)

    delta = (a0 + alpha) / 2
    if not delta.is_square():
        delta -= alpha

    x0 = new_sqrt_Fp(delta, exp=exp)
    x1 = a1 / (x0 + x0)

    return F([x0, x1])


def sqrt_Fp2(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4.
    """
    Fp2 = a.parent()
    p = Fp2.characteristic()
    i = Fp2.gen()  # i = √-1

    # Removing these asserts will speed things
    # up but adding them now to avoid annoying 
    # bugs
    assert p % 4 == 3
    assert i*i == -1

    a1 = a ** ((p - 3) // 4)
    x0 = a1 * a
    alpha = a1 * x0

    if alpha == -1:
        x0 *= i
    else:
        b = (1 + alpha) ** ((p - 1) // 2)
        x0 *= b

    return x0


def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)
