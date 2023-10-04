import time

from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL, ThetaCGLRadical4, ThetaCGLRadical8
from dim2 import ThetaCGLDim2
from dim3 import ThetaCGLDim3

from utilities import sqrt_Fp2, new_sqrt_Fp2, print_info

def time_function_ns(f):
    t0 = time.process_time_ns()
    eval(f)
    return (time.process_time_ns() - t0)

def time_ms(f):
    v = time_function_ns(f)
    return v / 1_000_000

def check(O0, O1, O2):
    print(f"Are the isogeneous curves isogeneous? {O0.cardinality() == O1.cardinality()}, {O0.cardinality() == O2.cardinality()}")
    print(f"Are the isogeneous curves the same? {O1.j_invariant() == O2.j_invariant()}, {O0.j_invariant(), O1.j_invariant(), O2.j_invariant()}")


p = 4 * 2**72 * 3**41 - 1  # any p = 3 mod 4
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
m1 = [1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1]
m2 = [0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1]

print_info(f"Example in dim 1")
print("- Sanity checks")
O0 = ThetaCGL(E0)
O1 = O0.bit_string(m1)
O2 = O0.bit_string(m2)
print(f"First message gives {O1} which hashes to {O1.to_hash()}")
print(f"Second message gives {O2} which hashes to {O2.to_hash()}")
check(O0, O1, O2)

print("- SageMath only")
O0 = ThetaCGL(E0)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print("- Faster sqrt function")
O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print("- New sqrt function")
O0 = ThetaCGL(E0, sqrt_function=new_sqrt_Fp2)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print_info(f"Example of a 4-radical isogeny")
O0 = ThetaCGLRadical4(E0, zeta4=F.gen())

print("- Sanity checks")
O1 = O0.bit_string(m1)
O2 = O0.bit_string(m2)
check(O0, O1, O2)

print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print_info(f"Example of a 8-radical isogeny")
i = F.gen()
zeta = i.sqrt()
sqrt2 = F(2).sqrt()
O0 = ThetaCGLRadical8(E0, zeta8=zeta, zeta4 = i, sqrt2 = sqrt2)

print("- Sanity checks")
O1 = O0.bit_string(m1)
O2 = O0.bit_string(m2)
check(O0, O1, O2)

print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print_info(f"Timings")
O0 = ThetaCGL(E0)
t_sage = time_ms("O0.hash(m1)")
print(f"Sage Sqrt Hashing time took: {t_sage}ms")

O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
t_fast = time_ms("O0.hash(m1)")
print(f"Fast Sqrt Hashing time took: {t_fast}ms")

O0 = ThetaCGL(E0, sqrt_function=new_sqrt_Fp2)
t_fast = time_ms("O0.hash(m1)")
print(f"New Sqrt Hashing time took: {t_fast}ms")

print_info(f"Example in dim 2")

O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")
m3 = [0, 1, 1, 0, 1, 0, 1] #length non multiple of 3
print(f"Hashing test 3: {O0.hash(m3)}")

O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=sqrt_Fp2)
print(f"Hashing test 3 (sqrt_Fp2): {O0.hash(m3)}")
O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=new_sqrt_Fp2)
print(f"Hashing test 3 (new_sqrt_Fp2): {O0.hash(m3)}")

print_info(f"Dimension Two Timings")
O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0)
t_sage = time_ms("O0.hash(m1)")
print(f"Sage Sqrt Hashing time took: {t_sage}ms")

O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=sqrt_Fp2)
t_fast = time_ms("O0.hash(m1)")
print(f"Fast Sqrt Hashing time took: {t_fast}ms")

O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0, sqrt_function=new_sqrt_Fp2)
t_fast = time_ms("O0.hash(m1)")
print(f"New Sqrt Hashing time took: {t_fast}ms")



print_info(f"Example in dim 3")
m1 = [1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1]
m2 = [0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1]
O0 = ThetaCGLDim3.from_elliptic_curves(E0, E0, E0)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")




def bench(N):
    for _ in range(N):
        O0.hash(m1)
    return
