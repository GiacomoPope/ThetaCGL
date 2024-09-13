"""
Symbolic verification of Proposition 20 in the paper 

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
(u0 : u1) is an 8-torsion point satisfying 
(a0^2 : a1^2) = (u0^4 + u1^4 : 2 * u0^2 * u1^2)

lam = (u0^8 - u1^8)
"""
relation1 = (u0^4 + u1^4) * a1^2 - 2 * u0^2 * u1^2 * a0^2
relation2 = lam^8 - (u0^8 - u1^8)
relation3 = sqrt2^2 - 2

I = K.ideal([relation1, relation2,relation3])

b0 = u0^2 + lam^2
b1 = u0^2 - lam^2

#v0 = sqrt(2) * a0 * (u0^2 - lam^2)
#v1 = a0^2 + lam^4 - 2*a0*u0*lam
v0 = a0*a1*(u0^2-lam^2)
v1 = a0^2*u0*u1 + lam^4*a1^2/(2*u0*u1) - sq2*lam*a0*a1*u0


# we prove that (v0 : v1) is an 8-torsion point of the correct form 
# i.e. (b0^2 : b1^2) = (v0^4 + v1^4 : 2 * v0^2 * v1^2)
relation = (v0^4 + v1^4) * b1^2 - 2 * v0^2 * v1^2 * b0^2

assert relation.numerator().reduce(I) == 0