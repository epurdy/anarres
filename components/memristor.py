from util import get_next_id
from components.component import TwoNodeCircuitComponent


class Memristor(TwoNodeCircuitComponent):
    def __init__(self, crkt, R_on, R_off, C=1e-3):
        self.id = get_next_id('M')
        super(Memristor, self).__init__(crkt)
        self.R_on = R_on
        self.R_off = R_off

        # a capacitor in parallel to model parasitic capcitance
        self.parasitic_capacitor = Capacitor(crkt, C)
        self.parasitic_capacitor.positive = self.positive
        self.parasitic_capacitor.negative = self.negative

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return True

    def M(self, x):
        return self.R_on * (1 - x) + self.R_off * x

    def Minv(self, x):
        return 1.0 / self.M(x)

    def negMinv(self, x):
        return -1.0 / self.M(x)

    def fM(self, pos, neg, h):
        i = (pos - neg) * self.Minv(h)
        return dcm.derivative_of_memristor_variables(pos - neg, h)

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp2(posvar, posvar, self.Minv,    [hvar]),
            FuncStamp2(posvar, negvar, self.negMinv, [hvar]),
            FuncStamp2(negvar, posvar, self.negMinv, [hvar]),
            FuncStamp2(negvar, negvar, self.Minv,    [hvar]),
        ]

    def Bstamp(self):
        hvar = Var(self.id, 'h')
        return [
            ConstStamp2(hvar, hvar, 1)
        ]

    def cstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp1(hvar, self.fM, [posvar, negvar, hvar])
        ]