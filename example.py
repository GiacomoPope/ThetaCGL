import time

from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL
from dim2 import ThetaCGLDim2
from utilities import sqrt_Fp2, new_sqrt_Fp2, print_info

def time_function_ns(f):
    t0 = time.process_time_ns()
    eval(f)
    return (time.process_time_ns() - t0)

def time_ms(f):
    v = time_function_ns(f)
    return v / 1_000_000

p = 4 * 2**72 * 3**41 - 1  # any p = 3 mod 4
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
m1 = [0, 1, 1, 0, 1, 1]
m2 = [0, 1, 1, 0, 1, 0]

print_info(f"Example")
print("SageMath only")
O0 = ThetaCGL(E0)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print("Faster sqrt function")
O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print("New sqrt function")
O0 = ThetaCGL(E0, sqrt_function=new_sqrt_Fp2)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print_info(f"Example in dim 2")
O0 = ThetaCGLDim2.from_elliptic_curves(E0, E0)
print(f"Hashing test 1: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")
m3 = [0, 1, 1, 0, 1, 0, 1] #length non multiple of 3
print(f"Hashing test 3: {O0.hash(m3)}")

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

def bench(N):
    for _ in range(N):
        O0.hash(m1)
    return
