#radical (4,4)-isogenies
def dim2_isogeny(a,b,c,d):
	a1 = a^2
	b1 = b^2
	c1 = c^2
	d1 = d^2

	AA = a1 + b1 + c1 + d1
	BB = a1 - b1 + c1 - d1
	CC = a1 + b1 - c1 - d1
	DD = a1 - b1 - c1 + d1

	newa = AA + sqrt(AA*BB) + sqrt(AA*CC) + sqrt(AA*DD)
	newb = AA - sqrt(AA*BB) + sqrt(AA*CC) - sqrt(AA*DD)
	newc = AA + sqrt(AA*BB) - sqrt(AA*CC) - sqrt(AA*DD)
	newd = AA - sqrt(AA*BB) - sqrt(AA*CC) + sqrt(AA*DD)

	return newa, newb,newc,newd

def chain(a,b,c,d, n=2):
	for i in range(n):
		a,b,c,d = dim2_isogeny(a,b,c,d)
	return a,b,c,d


[a,b,c,d] = var("a,b,c,d")

a4,b4,c4,d4 = chain(a,b,c,d)

new_a4 = 1
new_b4 = (b4/a4).canonicalize_radical()
new_c4 = (c4/a4).canonicalize_radical()
new_d4 = (d4/a4).canonicalize_radical()

den = new_b4.denominator()
new_a4 = new_a4*den
new_b4 = new_b4*den 
new_c4 = new_c4*den
new_d4 = new_d4*den

#simplify using the AA,BB,CC,DD 
AA = a^2 + b^2 + c^2 + d^2
BB = a^2 - b^2 + c^2 - d^2
CC = a^2 + b^2 - c^2 - d^2
DD = a^2 - b^2 - c^2 + d^2

new_a4 - (+sqrt(2)*sqrt(sqrt(AA)*sqrt(CC) + sqrt(BB)*sqrt(DD)) + sqrt(2)*sqrt(sqrt(AA)*sqrt(BB) + sqrt(CC)*sqrt(DD)) + sqrt(2)*sqrt(sqrt(AA)*sqrt(DD) + sqrt(BB)*sqrt(CC)) + 2*a)
new_b4 - (+sqrt(2)*sqrt(sqrt(AA)*sqrt(CC) + sqrt(BB)*sqrt(DD)) - sqrt(2)*sqrt(sqrt(AA)*sqrt(BB) + sqrt(CC)*sqrt(DD)) - sqrt(2)*sqrt(sqrt(AA)*sqrt(DD) + sqrt(BB)*sqrt(CC)) + 2*a)
new_c4 - (-sqrt(2)*sqrt(sqrt(AA)*sqrt(CC) + sqrt(BB)*sqrt(DD)) + sqrt(2)*sqrt(sqrt(AA)*sqrt(BB) + sqrt(CC)*sqrt(DD)) - sqrt(2)*sqrt(sqrt(AA)*sqrt(DD) + sqrt(BB)*sqrt(CC)) + 2*a)
new_d4 - (-sqrt(2)*sqrt(sqrt(AA)*sqrt(CC) + sqrt(BB)*sqrt(DD)) - sqrt(2)*sqrt(sqrt(AA)*sqrt(BB) + sqrt(CC)*sqrt(DD)) + sqrt(2)*sqrt(sqrt(AA)*sqrt(DD) + sqrt(BB)*sqrt(CC)) + 2*a)

########################################################

[A,B,C,D] = var("A,B,C,D")
term1 = sqrt(sqrt(A*B) + sqrt(C*D))
term2 = sqrt(sqrt(A*C) + sqrt(B*D))
term3 = sqrt(sqrt(A*D) + sqrt(B*C))

alpha1 = (2*sqrt(A*B*C*D) + A*B + C*D)^(1/4)
alpha2 = (2*sqrt(A*B*C*D) + A*C + B*D)^(1/4)
alpha3 = (2*sqrt(A*B*C*D) + C*B + A*D)^(1/4)

(term1^4 - alpha1^4).canonicalize_radical() == 0
(term2^4 - alpha2^4).canonicalize_radical() == 0
(term3^4 - alpha3^4).canonicalize_radical() == 0

