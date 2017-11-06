from anarres.util import get_next_id
from anarres.components.component import TwoNodeCircuitComponent
from anarres.components.vcvg import Vcvg
from anarres.components.memristor import Memristor
from anarres.components.resistor import Resistor


class DynamicCorrectionModule(TwoNodeCircuitComponent):
    """Dynamic correction module, as seen in Traversa and
    di Ventra (2015)"""
    def __init__(self, crkt, params, nodes, r_on, r_off):
        self.id = get_next_id('DCM')
        super(DynamicCorrectionModule, self).__init__(crkt)

        gnd = crkt._ground

        self.negative = gnd

        self.vcvg_b = Vcvg(crkt, params[0], nodes)
        self.vcvg_c = Vcvg(crkt, params[1], nodes)
        self.vcvg_d = Vcvg(crkt, params[4], nodes)
        self.vcvg_e = Vcvg(crkt, params[2], nodes)
        self.vcvg_f = Vcvg(crkt, params[3], nodes)

        self.mem_b = Memristor(crkt, r_on, r_off)
        self.mem_c = Memristor(crkt, r_on, r_off)
        self.mem_e = Memristor(crkt, r_on, r_off)
        self.mem_f = Memristor(crkt, r_on, r_off)

        self.res_d = Resistor(crkt, r_off)

        # hook negative of vcvg's to ground
        self.vcvg_b.negative = gnd
        self.vcvg_c.negative = gnd
        self.vcvg_d.negative = gnd
        self.vcvg_e.negative = gnd
        self.vcvg_f.negative = gnd

        # hook positive of vcvg's to negative of corresponding
        # memristor/resistor
        self.vcvg_b.positive = self.mem_b.negative
        self.vcvg_c.positive = self.mem_c.negative
        self.vcvg_d.positive = self.res_d.negative
        self.vcvg_e.positive = self.mem_e.negative
        self.vcvg_f.positive = self.mem_f.negative
        
        # hook positive of memristor/resistors to positive output
        self.mem_b.positive = self.positive
        self.mem_c.positive = self.positive
        self.res_d.positive = self.positive
        self.mem_e.positive = self.positive
        self.mem_f.positive = self.positive

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False
