from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL
from utilities import sqrt_Fp2

p=4*2**72*3**41-1 #any p = 3 mod 4
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
O0 = ThetaCGL(E0)
print(O0)
O1 = O0.advance()
print(O1)
print(O1.to_hash())

m1=[0,1,1,0,1,1]
m2=[0,1,1,0,1,0]
print(f"Hashing test: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")

print("Faster sqrt function")
O0 = ThetaCGL(E0, sqrt_function=sqrt_Fp2)
O1 = O0.advance()
print(O1)
print(O1.to_hash())
print(f"Hashing test: {O0.hash(m1)}")
print(f"Hashing test 2: {O0.hash(m2)}")
