import numpy as np

from anarres.util import get_next_id
from anarres.components.component import ThreeNodeCircuitComponent
from anarres.components.dcm import DynamicCorrectionModule
from anarres.components.resistor import Resistor

VC = 10
R_ON = 1.0
R_OFF = 2.0

class SoGate(ThreeNodeCircuitComponent):
    """
    v1--+   +--v2
        |   |
 G--DCM-+   +-DCM--G
        |   |
        R   R
        +_+_+
          |
  G--DCM--+
          |
          vo
    """
    def __init__(self, crkt, gate_params):
        self.id = get_next_id('DCM')
        super(SoGate, self).__init__(crkt)

        gnd = crkt._ground

        self.dcm1 = DynamicCorrectionModule(crkt, gate_params[0], [self.node1, self.node2, self.node3], 
                                            r_on=R_ON, r_off=R_OFF)
        self.dcm2 = DynamicCorrectionModule(crkt, gate_params[1], [self.node1, self.node2, self.node3],
                                            r_on=R_ON, r_off=R_OFF)
        self.dcmo = DynamicCorrectionModule(crkt, gate_params[2], [self.node1, self.node2, self.node3],
                                            r_on=R_ON, r_off=R_OFF)

        self.r1 = Resistor(crkt, R_OFF)
        self.r2 = Resistor(crkt, R_OFF)

        # hook up dcm's
        self.dcm1.positive = self.node1
        self.dcm2.positive = self.node2
        self.dcmo.positive = self.node3

        # hook up resistors
        self.r1.negative = self.node1
        self.r1.positive = self.node3
        self.r2.negative = self.node2
        self.r2.positive = self.node3

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

GATE_PARAMS = dict(
    AND=np.array(
        [[[ 0, -1,  1,      VC],  # B
          [ 1,  0,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 1,  0,  0,       0],  # F
          [ 4,  1, -3,     -VC]], # D

         [[-1,  0,  1,      VC],  # B
          [ 0,  1,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 0,  1,  0,       0],  # F
          [ 1,  4, -3,     -VC]], # D

         [[ 1,  0,  0,       0],  # B
          [ 0,  1,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 2,  2, -1, -2 * VC],  # F
          [-4, -4,  7,  2 * VC]], # D
     ]),

    OR=np.array(
        [[[ 0,  0,  1,       0],  # B
          [ 1,  0,  0,       0],  # C
          [ 0, -1,  1,     -VC],  # E
          [ 1,  0,  0,       0],  # F
          [ 4,  1, -3,      VC]], # D
         
         [[ 0,  0,  1,       0],  # B
          [ 0,  1,  0,       0],  # C
          [-1,  0,  1,     -VC],  # E
          [ 0,  1,  0,       0],  # F
          [ 1,  4, -3,      VC]], # D

         [[ 0,  0,  1,       0],  # B
          [ 2,  2, -1,  2 * VC],  # C
          [ 1,  0,  0,       0],  # E
          [ 0,  1,  0,       0],  # F
          [-4, -4,  7, -2 * VC]], # D
     ]),

    XOR=np.array(
        [[[ 0, -1, -1,      VC],  # B
          [ 0,  1,  1,      VC],  # C
          [ 0, -1,  1,     -VC],  # E
          [ 0,  1, -1,     -VC],  # F
          [ 6,  0, -1,       0]], # D

         [[-1,  0, -1,      VC],  # B
          [ 1,  0,  1,      VC],  # C
          [-1,  0,  1,     -VC],  # E
          [ 1,  0, -1,     -VC],  # F
          [ 0,  6, -1,       0]], # D

         [[-1, -1,  0,      VC],  # B
          [ 1,  1,  0,      VC],  # C
          [-1,  1,  0,     -VC],  # E
          [ 1, -1,  0,     -VC],  # F
          [-1, -1,  7,       0]], # D
     ]),
)

class SoAnd(SoGate):
    def __init__(self, crkt):
        super(SoAnd, self).__init__(crkt, GATE_PARAMS['AND'])

class SoOr(SoGate):
    def __init__(self, crkt):
        super(SoOr, self).__init__(crkt, GATE_PARAMS['OR'])

class SoXor(SoGate):
    def __init__(self, crkt):
        super(SoXor, self).__init__(crkt, GATE_PARAMS['XOR'])
