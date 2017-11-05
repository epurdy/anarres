import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import circuit2


def test_and_gate():
    crkt = circuit2.Circuit()
    and_gate = circuit2.SoAnd(crkt)

    v1_aux = circuit2.Vcdcg(crkt)
    v2_aux = circuit2.Vcdcg(crkt)
    v3_aux = circuit2.Vcdcg(crkt)

    and_gate.node1 = v1_aux.positive
    and_gate.node2 = v2_aux.positive
    and_gate.node3 = v3_aux.positive

    v1_aux.negative = v2_aux.negative = v3_aux.negative = crkt._ground

    # hardwire_v1 = circuit2.VoltageSource(crkt, -10)
    # hardwire_v1.positive = and_gate.node1
    # hardwire_v1.negative = crkt._ground

    crkt.normalize()
    print crkt.graph.edges()
    mna = crkt.assemble_mna_equation()
    stuff = mna.simulate_be(10.0, 0.1, 
                            vars=[circuit2.Var(and_gate.node1, 'v'),
                                  circuit2.Var(and_gate.node2, 'v'),
                                  circuit2.Var(and_gate.node3, 'v')])

    plt.plot(stuff)
    plt.savefig('and_dyn.png')

    # crkt.draw()
    # plt.savefig('and.png')

class Solc(object):
    def __init__(self):
        self.variables = dict()
        self.circuit = circuit2.Circuit()

    def add_variable(self, name):
        aux = circuit2.Vcdcg(self.circuit)
        aux.negative = self.circuit._ground
        self.variables[name] = aux.positive
        return name

    def and_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = circuit2.SoAnd(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def or_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = circuit2.SoOr(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def xor_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = circuit2.SoXor(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def set_on(self, name):
        hardwired_voltage = circuit2.VoltageSource(self.circuit, 10)        
        hardwired_voltage.positive = self.circuit._ground
        hardwired_voltage.negative = self.variables[name]

    def set_off(self, name):
        hardwired_voltage = circuit2.VoltageSource(self.circuit, -10)
        hardwired_voltage.positive = self.circuit._ground
        hardwired_voltage.negative = self.variables[name]

def test_factorization():
    solc = Solc()

    a0 = solc.add_variable('a0')
    a1 = solc.add_variable('a1')
    b0 = solc.add_variable('b0')
    b1 = solc.add_variable('b1')

    #        A1   A0
    #        B1   B0
    # --------------
    #      A1B0 A0B0
    # A1B1 A0B1

    ab0 = solc.and_gate('ab0', a0, b0)
    ab1_1 = solc.and_gate('ab1_1', a0, b1)
    ab1_2 = solc.and_gate('ab1_2', a1, b0)
    ab1 = solc.xor_gate('ab1', ab1_1, ab1_2)
    ab1_carry = solc.and_gate('ab1_carry', ab1_1, ab1_2)
    ab2_1 = solc.and_gate('ab2_1', a1, b1)
    ab2 = solc.xor_gate('ab2', ab2_1, ab1_carry)
    ab3 = solc.and_gate('ab3', ab2_1, ab1_carry)

    # set output to 4 = 0 1 0 0
    # set output to 9 = 1 0 0 1
    # set output to 7 = 0 1 1 1
    n = 4
    
    if n & 8:
        solc.set_on(ab3)
    else:
        solc.set_off(ab3)
    if n & 4:
        solc.set_on(ab2)
    else:
        solc.set_off(ab2)
    if n & 2:
        solc.set_on(ab1)
    else:
        solc.set_off(ab1)
    if n & 1:
        solc.set_on(ab0)
    else:
        solc.set_off(ab0)


    mna = solc.circuit.assemble_mna_equation()

    # x0 = np.concatenate(
    #     [2 * np.random.random(mna.nv) - 1,
    #      2 * np.random.random(mna.ni) - 1,
    #      np.random.random(mna.nh)])
    #x0 = None
    x0 = np.zeros(mna.n)

    stuff = mna.simulate_be(0.1, 0.001,
                            x0=x0,
                            vars=[circuit2.Var(solc.variables[a0], 'v'),
                                  circuit2.Var(solc.variables[a1], 'v'),
                                  circuit2.Var(solc.variables[b0], 'v'),
                                  circuit2.Var(solc.variables[b1], 'v'),
                              ])

    a0_val = stuff[-1, 0] >= 0
    a1_val = stuff[-1, 1] >= 0
    b0_val = stuff[-1, 2] >= 0
    b1_val = stuff[-1, 3] >= 0

    print 'A = ', 2 * a1_val + a0_val
    print 'B = ', 2 * b1_val + b0_val

    # raw_input()

    # stuff = mna.simulate_be(1.0, 0.001, 
    #                         x0=stuff[-1],
    #                         vars=[circuit2.Var(solc.variables[a0], 'v'),
    #                               circuit2.Var(solc.variables[a1], 'v'),
    #                               circuit2.Var(solc.variables[b0], 'v'),
    #                               circuit2.Var(solc.variables[b1], 'v'),
    #                               ])

    plt.plot(stuff[9 * len(stuff) // 10:])
    plt.savefig('factorization.png')

    # crkt.draw()
    # plt.savefig('and.png')

if __name__ == '__main__':
    test_factorization()
