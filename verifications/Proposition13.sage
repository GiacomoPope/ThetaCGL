"""
Radical 8-isogeny formula in dimension 1:

Symbolic verification of Proposition 13 in the paper 

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


K.<lam,u0,u1,a0,a1,sqrt2> = QQ[]
"""
(a0 : a1) is a theta nullpoint
(u0 : u1) is an 8-torsion point 
lam = (u0^8 - u1^8)
"""

relation1 = lam^8 - (u0^8 - u1^8)
relation2 = sqrt2^2 - 2

"""
The following relation describes the fact that
 2 * (u0 : u1)  = (1 : 0) 
"""

relation3 = (u0^4 + u1^4) * a1^2 - 2 * u0^2 * u1^2 * a0^2

I = K.ideal([relation1, relation2,relation3])

"""
we show that (b0 : b1) as below is the theta nullpoint of the codomain of an 8-isogeny
"""
b0 = u0^2 + lam^2
b1 = u0^2 - lam^2

v0 = a0*a1*(u0^2-lam^2)
v1 = a0^2*u0*u1 + lam^4*a1^2/(2*u0*u1) - sq2*lam*a0*a1*u0


# we prove that (v0 : v1) is an 8-torsion point of the correct form 
# i.e. (b0^2 : b1^2) = (v0^4 + v1^4 : 2 * v0^2 * v1^2)
relation = (v0^4 + v1^4) * b1^2 - 2 * v0^2 * v1^2 * b0^2

assert relation.numerator().reduce(I) == 0