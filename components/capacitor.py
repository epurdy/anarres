from util import get_next_id
from components.component import TwoNodeCircuitComponent
from circuit import Var, ConstStamp2, ConstStamp1


class Capacitor(TwoNodeCircuitComponent):
    """A capacitor without an explicit current variable.
    """
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor, self).__init__(crkt)
        self.capacitance = float(capacitance)

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

    def Bstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        return [
            ConstStamp2(posvar, posvar,  self.capacitance),
            ConstStamp2(posvar, negvar, -self.capacitance),
            ConstStamp2(negvar, posvar, -self.capacitance),
            ConstStamp2(negvar, negvar,  self.capacitance),
        ]


class Capacitor2(TwoNodeCircuitComponent):
    """A capacitor with an explicit current variable.
    """
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor2, self).__init__(crkt)
        self.capacitance = float(capacitance)

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self, static=False):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')

        if static:
            # capacitors become open circuits in static case
            return []
        else:
            return [
                # these do not appear in Najm but seem to be necessary
                ConstStamp2(posvar, ivar,  1),
                ConstStamp2(negvar, ivar, -1),

                # this one does appear in Najm
                ConstStamp2(ivar, ivar, 1),
            ]

    def Bstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(ivar, posvar, -self.capacitance),
            ConstStamp2(ivar, negvar, self.capacitance),
        ]
