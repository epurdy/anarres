import numpy as np

from anarres.util import get_next_id
from anarres.components.component import TwoNodeCircuitComponent
from anarres.circuit import Var, ConstStamp2, FuncStamp1
import anarres.vcdcg as vcdcg

class Vcdcg(TwoNodeCircuitComponent):
    """Voltage-controlled differential current generator, as seen in Traversa and
    di Ventra (2015)"""

    # since the Vcdcg's have collective behavior, they all need to know about
    # one another
    _all_vcdcg_currents = []

    @classmethod
    def all_vcdcg_currents(cls):
        return cls._all_vcdcg_currents

    def __init__(self, crkt):
        self.id = get_next_id('VCDCG')
        super(Vcdcg, self).__init__(crkt)
        self.__class__._all_vcdcg_currents.append(Var(self.id, 'i'))

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return True

    def deriv_i(self, posv, negv, i, h):
        v = negv - posv
        rv = vcdcg.derivative_of_dcg_currents(v, i, h, delta_s=1e-3, gamma=1.0, v_c=10.0)
        return rv

    def deriv_h(self, h, *ivals):
        # note that the ivals are out of order but that it does not matter
        ivals = np.array(ivals)
        rv = vcdcg.derivative_of_dcg_variables(ivals, h, i_min=9.0, i_max=11.0, k_s=1.0, k_i=3.0, delta_i=1e-3)
        return rv

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
        ]

    def Bstamp(self):
        ivar = Var(self.id, 'i')
        hvar = Var(self.id, 'h')
        return [
            ConstStamp2(ivar, ivar, 1),
            ConstStamp2(hvar, hvar, 1),
        ]

    def cstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp1(ivar, self.deriv_i, [posvar, negvar, ivar, hvar]),
            FuncStamp1(hvar, self.deriv_h, [hvar] + self.all_vcdcg_currents())
        ]
