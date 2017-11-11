import anarres.solc as solc

class SoSat(object):
    def __init__(self, n):
        self.n = n

        self.circ = solc.Solc()

        self.vars = [self.circ.add_variable('a%d' % (i + 1)) for i in xrange(n)]

        self.circ.add_variable('true')
        self.circ.set_on('true')

    def add_literal(self, a):
        assert a != 0

        if a > 0:
            a = 'a%d' % a
        else:
            a = self.circ.xor_gate('a%d' % (-a), 'true')

        return a

    def add_clause(self, a, b, c):
        a = self.add_literal(a)
        b = self.add_literal(b)
        c = self.add_literal(c)
        
        aorb = self.circ.or_gate(a, b)
        aorborc = self.circ.or_gate(aorb, c)

        self.circ.set_on(aorborc)
