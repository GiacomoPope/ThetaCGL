from sage.all import GF, EllipticCurve
from dim1 import ThetaCGL

p=4*2**72*3**41-1
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
O0 = ThetaCGL(E0)
print(O0)
O1 = O0.advance()
print(O1)

print(f"Hashing test: {O0.hash([0,1,1,0])}")
