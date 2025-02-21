"""
Radical 4-isogeny formula in dimension 1:

Symbolic verification of Proposition 12 in the paper 

<< Radical 2-isogenies in theta coordinates and 
applications to cryptographic hash functions in
dimensions 1, 2 and 3 >>
"""

def square(c0,c1):
    return c0^2, c1^2

def hadamard(c0,c1):
    d0 = c0 + c1
    d1 = c0 - c1
    return d0, d1

def scale(c0,c1,d0,d1):
    return c0/d0, c1/d1

def radical_2_isogeny(a0,a1):
    x0, x1 = hadamard(a0 * a0, a1 * a1)
    x01 = x0 * x1
    y1 = sqrt(x01)
    b0, b1 = hadamard(x0, y1)
    return b0, b1

def radical_4_isogeny(a0,a1):
    x0, x1 = hadamard(a0 * a0, a1 * a1)
    x01 = x0 * x1
    lam = sqrt(sqrt(x01))
    b0,b1 = hadamard(a0, lam)
    return b0, b1

"""
We show that applying the 2-isogeny formula twice yields the same output 
as applying the radical 4-isogeny formula (assuming that the  choices of 
square / fourth roots are compatible).
"""

K.<a0,a1> = QQ[]

#radical 4-isogeny
b0,b1 = radical_4_isogeny(a0,a1)

# 2 radical 2-isogenies
B0,B1 = a0,a1
for i in range(2):
    B0,B1 = radical_2_isogeny(B0,B1)

(b0*B1 - b1*B0).canonicalize_radical() == 0