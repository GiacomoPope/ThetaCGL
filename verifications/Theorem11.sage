"""
Symbolic verification of Theorem 11 in the paper 

<< Radical 2-isogenies in theta coordinates and 
applications to cryptographic hash functions in
dimensions 1, 2 and 3 >>

"""


#We consider an isogeny phi: A -> B, and use the same notation as in the proof:
# a = (a000, ..., a111) the theta null point on A
# x = (x000, ..., x111) = Hadamard ( Square ( a ) )
# y = (y000, ..., y111) = Ts(x), i.e. the dual theta null point of B

# The goal is to compute y111 given y000, ... , y100 and a000, ..., a111

K.<y000,y001,y010,y011,y100,y101,y110,y111> = QQ[]
# We note that y is a theta null point, so satisfies the relation described in Proposition 10
# We compute beta, gamma, delta in terms of y:
beta00 = y000^2 + y001^2 + y010^2 + y011^2 + y100^2 + y101^2 + y110^2 + y111^2
beta01 = y000^2 - y001^2 + y010^2 - y011^2 + y100^2 - y101^2 + y110^2 - y111^2
beta10 = y000^2 + y001^2 - y010^2 - y011^2 + y100^2 + y101^2 - y110^2 - y111^2
beta11 = y000^2 - y001^2 - y010^2 + y011^2 + y100^2 - y101^2 - y110^2 + y111^2

delta00 = y000^2 + y001^2 + y010^2 + y011^2 - y100^2 - y101^2 - y110^2 - y111^2
delta01 = y000^2 - y001^2 + y010^2 - y011^2 - y100^2 + y101^2 - y110^2 + y111^2
delta10 = y000^2 + y001^2 - y010^2 - y011^2 - y100^2 - y101^2 + y110^2 + y111^2
delta11 = y000^2 - y001^2 - y010^2 + y011^2 - y100^2 + y101^2 + y110^2 - y111^2

gamma00 = 2*(y000*y100 + y001*y101 + y010*y110 + y011*y111)
gamma01 = 2*(y000*y100 - y001*y101 + y010*y110 - y011*y111)
gamma10 = 2*(y000*y100 + y001*y101 - y010*y110 - y011*y111)
gamma11 = 2*(y000*y100 - y001*y101 - y010*y110 + y011*y111)

#We compute R1, R2, R3 and express the relation from Proposition 10.
R1 = beta00 * beta01 * beta10 * beta11
R2 = gamma00 * gamma01 * gamma10 * gamma11
R3 = delta00 * delta01 * delta10 * delta11

relation = R1^2 + R2^2 + R3^2 - 2*(R1*R2 + R1*R3 + R2*R3)
I = K.ideal(relation)

# Now we want to express things in terms of x000, ..., x111 as in Theorem 11
# Note that yi^2 = xi*x000 for all i
x000 = y000
[x001, x010, x011, x100, x101, x110, x111] = [y001^2/x000, y010^2/x000, y011^2/x000, y100^2/x000, y101^2/x000, y110^2/x000, y111^2/x000]

# T from Theorem 11
T = 16*(x000*x100 - x001*x101 + x010*x110 - x011*x111)^2 - 64*(x000*x010*x100*x110 + x001*x011*x101*x111)

#we note that the beta_i and gamma_i in the theorem are also scaled by x000:
#that is: betai_theorem = betai * x000 for all i, and in particular:
R1_theorem = R1 / x000^4
R3_theorem = R3 / x000^4


#######################################
##### CASE 1: R1 + R3 - T nonzero #####
#######################################

# Claim 1 in the proof:
assert T == (R2 - 2^7 * y000*y001*y010*y011*y100*y101*y110*y111)/x000^4

# formula for y111 = y_numerator / y_denominator
y_numerator = x000^3 * ((R1_theorem + R3_theorem - T)^2 + 2^14 * x000*x001*x010*x011*x100*x101*x110*x111 - 4*R1_theorem*R3_theorem)
y_denominator = 2^8 * (R1_theorem + R3_theorem - T) * y001*y010*y011*y100*y101*y110

assert (y_numerator - y111*y_denominator).numerator().reduce(I) == 0

####################################
##### CASE 2: R1 + R3 - T = 0  #####
####################################

relation2 = (R1_theorem + R3_theorem - T).numerator()
J = I + K.ideal(relation2)

# Here, we only show that the formula is correct up to sign. 
# For the choice of sign, we refer to the proof in the paper.

# We denote by A000, ... , A111 the squares of the theta null coordinates of A
A000 = 1/2^3 * (x000 + x001 + x010 + x011 + x100 + x101 + x110 + x111)
A001 = 1/2^3 * (x000 - x001 + x010 - x011 + x100 - x101 + x110 - x111)
A010 = 1/2^3 * (x000 + x001 - x010 - x011 + x100 + x101 - x110 - x111)
A011 = 1/2^3 * (x000 - x001 - x010 + x011 + x100 - x101 - x110 + x111)
A100 = 1/2^3 * (x000 + x001 + x010 + x011 - x100 - x101 - x110 - x111)
A101 = 1/2^3 * (x000 - x001 + x010 - x011 - x100 + x101 - x110 + x111)
A110 = 1/2^3 * (x000 + x001 - x010 - x011 - x100 - x101 + x110 + x111)
A111 = 1/2^3 * (x000 - x001 - x010 + x011 - x100 + x101 + x110 - x111)

y_numerator = (-2^6 * x000^3)^2 * A000*A001*A010*A011*A100*A101*A110*A111
y_denominator = (y001*y010*y011*y100*y101*y110)^2

assert (y_numerator - y111^2*y_denominator).numerator().reduce(J) == 0

# we further verify some computational claims in the proof. 

# claim 2
assert (2^14 * x000*x001*x010*x011*x100*x101*x110*x111 - 4 * R1_theorem * R3_theorem).numerator().reduce(J) == 0

# claim 3
assert (x000*x001*x010*x011*x100*x101*x110*x111 - 2^12*A000*A001*A010*A011*A100*A101*A110*A111).numerator().reduce(J) == 0