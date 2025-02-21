"""
Radical 4-isogeny formula in dimension 2:

Symbolic verification of Proposition 14 in the paper 

<< Radical 2-isogenies in theta coordinates and 
applications to cryptographic hash functions in
dimensions 1, 2 and 3 >>
"""

def hadamard(x,y,z,t):
    return (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t)

def radical_2_isogeny(a0,a1,a2,a3):
	x0,x1,x2,x3 = hadamard(a0 * a0, a1 * a1, a2 * a2, a3 * a3)
	b0,b1,b2,b3 = hadamard(x0, sqrt(x0)*sqrt(x1), sqrt(x0)*sqrt(x2), sqrt(x0)*sqrt(x3))
	return b0,b1,b2,b3

a0,a1,a2,a3 = var("a0,a1,a2,a3")

# We apply the radical 2-isogeny formula 2 times.
B0,B1,B2,B3 = a0,a1,a2,a3
for i in range(2):
    B0,B1,B2,B3 = radical_2_isogeny(B0,B1,B2,B3)

# verification of the radical 4-isogeny divided into steps
x0,x1,x2,x3 = hadamard(a0 * a0, a1 * a1, a2 * a2, a3 * a3)
t1 = sqrt(sqrt(x0)*sqrt(x1) + sqrt(x2)*sqrt(x3))
t2 = sqrt(sqrt(x0)*sqrt(x2) + sqrt(x1)*sqrt(x3))
t3 = sqrt(sqrt(x0)*sqrt(x3) + sqrt(x1)*sqrt(x2))
bb0, bb1, bb2, bb3 = hadamard(sqrt(2)*a0, t1, t2, t3)

"""
we show that (bb0 : bb1 : bb2 : bb3) = (B0 : B1 : B2 : B3) 
"""
(bb0*B1 - bb1*B0).canonicalize_radical() == 0
(bb0*B2 - bb2*B0).canonicalize_radical() == 0
(bb0*B3 - bb3*B0).canonicalize_radical() == 0

"""
we use x0,x1,x2,x3 as variables 
having in mind x0,x1,x2,x3 = hadamard(a0 * a0, a1 * a1, a2 * a2, a3 * a3)

and we show that t1,t2,t2 = alpha1,alpha2,alpha3 as in the statement of
the proposition.
"""
x0, x1, x2, x3 = var("x0,x1,x2,x3")
t1 = sqrt(sqrt(x0)*sqrt(x1) + sqrt(x2)*sqrt(x3))
t2 = sqrt(sqrt(x0)*sqrt(x2) + sqrt(x1)*sqrt(x3))
t3 = sqrt(sqrt(x0)*sqrt(x3) + sqrt(x1)*sqrt(x2))

#bb0, bb1, bb2, bb3 = hadamard(sqrt(2)*a0, t1, t2, t3)

alpha14 = (2*sqrt(x0*x1*x2*x3) + x0*x1 + x2*x3)
alpha24 = (2*sqrt(x0*x1*x2*x3) + x0*x2 + x1*x3)
alpha34 = (2*sqrt(x0*x1*x2*x3) + x1*x2 + x0*x3)

(t1^4 - alpha14).canonicalize_radical() == 0
(t2^4 - alpha24).canonicalize_radical() == 0
(t3^4 - alpha34).canonicalize_radical() == 0