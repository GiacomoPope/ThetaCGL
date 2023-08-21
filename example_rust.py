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


p = 79*2**247 - 1
F = GF(p**2, name="i", modulus=[1, 0, 1])
# m1 = [1, 1, 1, 0, 1, 1, 1]
m1 = [1]



print_info(f"Example in dim 1")

null_coords = F([1,2]), F([2,1])
O0 = ThetaCGL(null_coords, sqrt_function=new_sqrt_Fp2)
out = O0.hash(m1)
print(out)
