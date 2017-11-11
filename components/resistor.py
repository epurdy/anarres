from util import get_next_id
from components.component import TwoNodeCircuitComponent
from circuit import Var, ConstStamp2, ConstStamp1


class Resistor(TwoNodeCircuitComponent):
    """A resistor without an explicit current variable.
    """
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor, self).__init__(crkt)
        self.resistance = float(resistance)

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

    def Astamp(self, static=False):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        return [
            ConstStamp2(posvar, posvar,  1.0 / self.resistance),
            ConstStamp2(posvar, negvar, -1.0 / self.resistance),
            ConstStamp2(negvar, posvar, -1.0 / self.resistance),
            ConstStamp2(negvar, negvar,  1.0 / self.resistance),
        ]


class Resistor2(TwoNodeCircuitComponent):
    """A resistor with an explicit current variable.
    """
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor2, self).__init__(crkt)
        self.resistance = float(resistance)

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self, static=False):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),
            ConstStamp2(ivar, ivar, -self.resistance)
        ]
