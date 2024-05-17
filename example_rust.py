import time

from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL, ThetaCGLRadical4
from dim2 import ThetaCGLDim2, ThetaCGLDim2Radical4
from dim3 import ThetaCGLDim3
from utilities import sqrt_Fp2, fourth_Fp2, print_info


def time_function_ns(f):
    t0 = time.process_time_ns()
    eval(f)
    return time.process_time_ns() - t0


def time_ms(f):
    v = time_function_ns(f)
    return v / 1_000_000


def to_hex_str(a):
    p = a.parent().characteristic()
    byte_len = (p.nbits() + 7) // 8
    a0, a1 = a.list()
    a0_bytes = int(a0).to_bytes(byte_len, byteorder="little")
    a1_bytes = int(a1).to_bytes(byte_len, byteorder="little")

    return (a0_bytes + a1_bytes).hex()


def to_little_u64(x):
    y = int(x)
    res = []
    while y:
        tmp = y % 2**64
        res.append(hex(tmp))
        y >>= 64
    l = fmt_little_u64(res)
    return fmt_list(l)


def to_little_u64_mont(x, n):
    x = 2 ** (64 * n) * x
    y = int(x)
    res = []
    while y:
        tmp = y % 2**64
        res.append(hex(tmp))
        y >>= 64
    l = fmt_little_u64(res)
    return fmt_list(l)


def fmt_little_u64(res):
    out = []
    for r in res:
        num = r[2:]
        new_num = num.zfill(16)
        new_num = new_num.upper()
        new_num = "0x" + new_num
        out.append(new_num)
    return out


def fmt_list(l):
    l = str(l)
    l = l.replace("[", "")
    l = l.replace("]", "")
    l = l.replace("'", "")
    return l


# fmt: off
m1 = [
    1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 
    0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 
    0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 
    0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 
    0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 
    1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1
]

if True:
    p = 79 * 2**247 - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E0 = EllipticCurve(F, [1, 0])


    print_info(f"Example dimension 1 radical 2 in p254")

    O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
    out = O0.hash(m1)
    print(out)

    a, b = O0.domain
    # print(f"Theta coordinates as hex strings: ")
    # print(to_hex_str(a))
    # print(to_hex_str(b))

    # print(f"Theta coordinates as u64 arrays (MONTGOMERY FORM): ")
    # print(to_little_u64_mont(a[0], 4))
    # print(to_little_u64_mont(a[1], 4))
    # print(to_little_u64_mont(b[0], 4))
    # print(to_little_u64_mont(b[1], 4))
    # print()

    print_info(f"Example dimension 1 radical 4 in p254")
    zeta = F.gen()
    O0 = ThetaCGLRadical4(E0, zeta4=zeta, fourth_root_function=fourth_Fp2)
    out = O0.hash(m1)
    print(out)

    p = 2**255 - 921
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    E0 = EllipticCurve(F, [1, 0])


    print_info(f"Example dimension 1 radical 2 in p921")
    O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
    out = O0.hash(m1)
    print(out)
    # a, b = O0.domain

    # print(f"Theta coordinates as hex strings: ")
    # print(to_hex_str(a))
    # print(to_hex_str(b))

    # print(f"Theta coordinates as u64 arrays: ")
    # print(to_little_u64(a[0]))
    # print(to_little_u64(a[1]))
    # print(to_little_u64(b[0]))
    # print(to_little_u64(b[1]))

    print_info(f"Example dimension 1 radical 4 in p921")
    zeta = F.gen()
    O0 = ThetaCGLRadical4(E0, zeta4=zeta, fourth_root_function=fourth_Fp2)
    out = O0.hash(m1)
    print(out)

if True:
    p = 2**127 - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    zeta = F.gen()
    E0 = EllipticCurve(F, [1, 0])

    print_info(f"Example dimension 2 radical 2 in p127")
    O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=sqrt_Fp2)
    out = O0.hash(m1)
    for h in out:
        print(h)

    # print(f"Theta coordinates as hex strings: ")
    # a, b, c, d = O0.domain
    # print(to_hex_str(a))
    # print(to_hex_str(b))
    # print(to_hex_str(c))
    # print(to_hex_str(d))

    # print(f"Theta coordinates as u64 arrays (MONTGOMERY FORM): ")
    # print(to_little_u64_mont(a[0], 2))
    # print(to_little_u64_mont(a[1], 2))
    # print(to_little_u64_mont(b[0], 2))
    # print(to_little_u64_mont(b[1], 2))
    # print(to_little_u64_mont(c[0], 2))
    # print(to_little_u64_mont(c[1], 2))
    # print(to_little_u64_mont(d[0], 2))
    # print(to_little_u64_mont(d[1], 2))

    print_info(f"Example dimension 2 radical 4 in p127")

    O0 = ThetaCGLDim2Radical4.from_elliptic_curves(E0, E0, sqrt_function=sqrt_Fp2, fourth_root_function=fourth_Fp2, zeta4=zeta)
    out = O0.hash(m1)
    for h in out:
        print(h)

if True:
    p = 2**64 - 257
    F = GF(p**2, name="i", modulus=[1, 0, 1])
    zeta = F.gen()
    E0 = EllipticCurve(F, [1, 0])

    print_info(f"Example dimension 3 radical 2 in p64")
    O0 = ThetaCGLDim3.from_elliptic_curves(E0, E0, E0)
    """
    out = O0.hash(m1)
    for h in out:
        print(h)
    """

    # print(f"Theta coordinates as hex strings: ")
    a, b, c, d, e, f, g, h = O0.domain
    print(to_hex_str(a))
    print(to_hex_str(b))
    print(to_hex_str(c))
    print(to_hex_str(d))
    print(to_hex_str(e))
    print(to_hex_str(f))
    print(to_hex_str(g))
    print(to_hex_str(h))

    print(f"Theta coordinates as u64 arrays (MONTGOMERY FORM): ")
    print(to_little_u64(a[0]))
    print(to_little_u64(a[1]))

    print(to_little_u64(b[0]))
    print(to_little_u64(b[1]))

    print(to_little_u64(c[0]))
    print(to_little_u64(c[1]))

    print(to_little_u64(d[0]))
    print(to_little_u64(d[1]))

    print(to_little_u64(e[0]))
    print(to_little_u64(e[1]))

    print(to_little_u64(f[0]))
    print(to_little_u64(f[1]))

    print(to_little_u64(g[0]))
    print(to_little_u64(g[1]))

    print(to_little_u64(h[0]))
    print(to_little_u64(h[1]))
