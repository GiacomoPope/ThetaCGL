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
        x = i * x0
    else:
        b = (1 + alpha) ** ((p - 1) // 2)
        x = b * x0

    return x


def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)
