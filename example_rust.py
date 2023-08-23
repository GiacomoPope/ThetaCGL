import time

from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL, ThetaCGLRadical4
from dim2 import ThetaCGLDim2
from utilities import sqrt_Fp2, new_sqrt_Fp2, print_info

def time_function_ns(f):
    t0 = time.process_time_ns()
    eval(f)
    return (time.process_time_ns() - t0)

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


p = 79*2**247 - 1
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
O0 = ThetaCGL(E0, sqrt_function=new_sqrt_Fp2)

"""
print_info(f"Example dimension 1 radical 2 in p254")

a, b = O0.domain
print(to_hex_str(a))
print(to_hex_str(b))

out = O0.hash(m1)
print(out)
"""

print_info(f"Example dimension 1 radical 4 in p254")
zeta = F.gen()
print("zeta:")
print(zeta)
O0 = ThetaCGLRadical4(E0, zeta=zeta)
out = O0.hash(m1)
print(out)

"""
print_info(f"Example in p127")
# p = 2**127 - 1
p = 27*2**122 - 1

F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=new_sqrt_Fp2)

a, b, c, d = O0.domain
print(to_hex_str(a))
print(to_hex_str(b))
print(to_hex_str(c))
print(to_hex_str(d))

out = O0.hash(m1)
for h in out:
    print(h)



# O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=new_sqrt_Fp2)

# a,b,c,d = O0.domain
# print(to_hex_str(a))
# print(to_hex_str(b))
# print(to_hex_str(c))
# print(to_hex_str(d))


# h = O0.hash(m1)
# print("hash:")
# print(h)
"""