"""
Factorize numbers using a self-organizing binary multiplier.
"""

import anarres.solc as solc


def so_adder(circ, a, b):
    n = max(len(a), len(b))

    a = a + [None] * (n - len(a))
    b = b + [None] * (n - len(b))
        
    answer = []
    carry = None
    for abit, bbit in zip(a, b):
        if abit is None and bbit is None and carry is None:
            assert False, 'should not happen'
        elif abit is None and carry is None:
            answer.append(bbit)
            carry = None
        elif bbit is None and carry is None:
            answer.append(abit)
            carry = None
        elif abit is None and bbit is None:
            answer.append(carry)
            carry = None
        elif abit is None:
            answer.append(circ.xor_gate(bbit, carry))
            carry = circ.and_gate(bbit, carry)
        elif bbit is None:
            answer.append(circ.xor_gate(abit, carry))
            carry = circ.and_gate(abit, carry)
        elif carry is None:
            answer.append(circ.xor_gate(abit, bbit))
            carry = circ.and_gate(abit, bbit)
        else:
            # a + b = x2 x1
            x1 = circ.xor_gate(abit, bbit)
            x2 = circ.and_gate(abit, bbit)

            # x2 x1 + carry = y2 y1
            y1 = circ.xor_gate(x1, carry)
            z = circ.and_gate(x1, carry)
            y2 = circ.or_gate(x2, z)

            carry = y2

    return answer + [carry]
                
            

class SoMultiplier(object):
    def __init__(self, n):
        self.n = n

        self.circ = solc.Solc()

        #             A2   A1   A0
        #             B2   B1   B0
        # ------------------------
        #           A2B0 A1B0 A0B0
        #      A2B1 A1B1 A0B1
        # A2B2 A2B1 A2B0

        self.a = []
        self.b = []
        for i in xrange(n):
            self.a.append(self.circ.add_variable('a%d' % i))
            self.b.append(self.circ.add_variable('b%d' % i))

        products = []
        for i in xrange(n):
            prod = [None] * i
            for j in xrange(n):
                ajbi = self.circ.and_gate(self.a[j], self.b[i])
                prod.append(ajbi)
            products.append(prod)

        print products

        total = products[0]
        for prod in products[1:]:
            print total, prod
            total = so_adder(self.circ, total, prod)

        self.answer = total

    def get_input_vars(self):
        return self.a, self.b

    def set_answer(self, N):
        self.N = N
        bits = list(bin(N)[2:])[::-1]
        bits = bits + [0] * (len(self.answer) - len(bits))

        if len(bits) > len(self.answer):
            assert False, 'N is too large (max=%d)' % (2 ** len(self.answer) - 1)

        for i, b in enumerate(bits):
            if b:
                self.circ.set_on(self.answer[i])
            else:
                self.circ.set_off(self.answer[i])
