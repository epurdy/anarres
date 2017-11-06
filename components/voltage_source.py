from util import get_next_id
from componets.component import TwoNodeCircuitComponent

    
class VoltageSource(TwoNodeCircuitComponent):
    """Voltage sources are required to have an explicit current variable."""
    def __init__(self, crkt, voltage):
        self.id = get_next_id('V')
        super(VoltageSource, self).__init__(crkt)
        self.voltage = voltage

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),
        ]

    def cstamp(self):
        ivar = Var(self.id, 'i')
        return [
            ConstStamp1(ivar, self.voltage),
        ]
