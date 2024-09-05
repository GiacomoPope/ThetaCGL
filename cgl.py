from utilities import sqrt_Fp2, fourth_Fp2


class CGL:
    def __init__(self, domain, chunk=1):
        self.domain = domain
        self.chunk = chunk

    def __repr__(self):
        return f"CGL with domain={self.domain}"

    def sqrt(self, x):
        return sqrt_Fp2(x)

    def fourth_root(self, x):
        return fourth_Fp2(x)

    def advance(self, bits=None):
        pass

    def bit_string(self, message):
        r = self
        for x in range(0, len(message), self.chunk):
            bits = message[x : x + self.chunk]
            if len(bits) < self.chunk:
                # pad bits if too short
                bits += [0] * (self.chunk - len(bits))
            r = r.advance(bits)
        return r

    def to_hash():
        pass

    def hash(self, bits):
        r = self.bit_string(bits)
        return r.to_hash()
