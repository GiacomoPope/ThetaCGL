# Defining equation for im(Th) subset PP^7

#Let (a0, ..., a7) be a level-2 theta nullpoint, then relation F holds (see Proposition 7 in our paper)
K.<a0,a1,a2,a3,a4,a5,a6,a7,a8> = QQ[]

#square
a02 = a0**2
a12 = a1**2
a22 = a2**2
a32 = a3**2
a42 = a4**2
a52 = a5**2
a62 = a6**2
a72 = a7**2

#hadamard
c0 = a02 + a12 + a22 + a32 + a42 + a52 + a62 + a72
c1 = a02 - a12 + a22 - a32 + a42 - a52 + a62 - a72
c2 = a02 + a12 - a22 - a32 + a42 + a52 - a62 - a72
c3 = a02 - a12 - a22 + a32 + a42 - a52 - a62 + a72
c4 = a02 + a12 + a22 + a32 - a42 - a52 - a62 - a72
c5 = a02 - a12 + a22 - a32 - a42 + a52 - a62 + a72
c6 = a02 + a12 - a22 - a32 - a42 - a52 + a62 + a72
c7 = a02 - a12 - a22 + a32 - a42 + a52 + a62 - a72

#remaining level-4 squares
b0 = 2*(a0*a4 + a1*a5 + a2*a6 + a3*a7)
b1 = 2*(a0*a4 - a1*a5 + a2*a6 - a3*a7)
b2 = 2*(a0*a4 + a1*a5 - a2*a6 - a3*a7)
b3 = 2*(a0*a4 - a1*a5 - a2*a6 + a3*a7)

R1 = c0*c1*c2*c3
R2 = b0*b1*b2*b3
R3 = c4*c5*c6*c7

F = R1**2 +  R2**2 + R3**2 - 2 * (R1*R2 + R1*R3 + R2*R3)


# Proof of Proposition 8
# for the first part, we express a7, in terms of a0,a1,a2,a3,a4,a5,a6  and a7^2 = a72

T = 16* ((a02*a42 - a12*a52 + a22*a62 - a32*a72)^2 - 4*(a02*a22*a42*a62 + a12*a32*a52*a72))
A7_num = (R1 + R3 - T)^2 + 2^14 * (a02*a12*a22*a32*a42*a52*a62*a72) - 4*R1*R3
A7_den =  2^8*(R1 + R3 - T)*a0*a1*a2*a3*a4*a5*a6

assert F == A7_num - a7*A7_den
#this provides us with a formula for a7:
A7 = A7_num / A7_den

#Note that in the situation of Proposition 8, we compute
#y7 = A7(y0,y1,y2,y3,y4,y5,y6,x7*x0).

# the above is only valid if the denominator is nonzero. 
# a0*a1*a2*a3*a4*a5*a6*a7 is nonzero (by assumption)
# now we prove the formula when R1 + R3 - T = 0

# this condition implies 2^12 * (a02*a12*a22*a32*a42*a52*a62*a72) - c0*c1*c2*c3*c4*c5*c6*c7 = 0:

condition1 = R1 + R3 - T
condition2 = 2^12 * (a02*a12*a22*a32*a42*a52*a62*a72) - c0*c1*c2*c3*c4*c5*c6*c7

assert condition2 in K.ideal(F, condition1)

#this means we can compute sqrt(c7) in terms of a0,a1,a2,a3,a4,a4,a6,a7