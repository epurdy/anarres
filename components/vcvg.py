from anarres.util import get_next_id
from anarres.components.component import TwoNodeCircuitComponent
from anarres.circuit import Var, ConstStamp2, FuncStamp1


class Vcvg(TwoNodeCircuitComponent):
    """Voltage-controlled voltage generator"""
    def __init__(self, crkt, params, nodes):
        self.id = get_next_id('VCVG')
        super(Vcvg, self).__init__(crkt)

        self.params = params
        self.nodes = nodes

    def normalize(self): 
        redo = super(Vcvg, self).normalize()       
        norm_nodes = [self._circuit.node_lut[node] for node in self.nodes]
        if norm_nodes != self.nodes:
            redo = True
            self.nodes = norm_nodes
        return redo

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def voltage_fn(self, *x):
        return self.params.dot(x + (1,))

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
        input_vars = [Var(node, 'v') for node in self.nodes]
        return [
            FuncStamp1(ivar, self.voltage_fn, input_vars),
        ]
