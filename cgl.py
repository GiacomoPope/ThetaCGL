from utilities import canonical_root


class CGL:
    def __init__(self, domain, sqrt_function=None, fourth_root_function=None, chunk=1):
        self.domain = domain
        self.sqrt_function = sqrt_function
        self.chunk = chunk

    def __repr__(self):
        return f"CGL with domain={self.domain}"

    def sqrt(self, x):
        if self.sqrt_function is None:
            r = x.sqrt()
        else:
            r = self.sqrt_function(x)

        return canonical_root(r)

    def fourth_root(self, x):
        if self.fourth_root_function is None:
            r = self.sqrt(self.sqrt(x))
        else:
            r = self.fourth_root_function(x)

        return canonical_root(r)

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
