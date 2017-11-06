import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from anarres.circuit import Circuit, Var
from anarres.components.all import Vcdcg, Resistor
import anarres.solc as solc

def test_vcdcgs():
    crkt = Circuit()
    vars = []
    for i in xrange(10):
        vcdcg = Vcdcg(crkt)
        res = Resistor(crkt, 1.0)
        vcdcg.negative = crkt._ground
        vcdcg.positive = res.negative
        res.positive = crkt._ground
        vars.append(Var(vcdcg.positive, 'v'))
        vars.append(Var(vcdcg.id, 'i'))
        vars.append(Var(vcdcg.id, 'h'))

    mna = crkt.assemble_mna_equation()
    x0 = np.concatenate([
        20 * np.random.random(mna.nv) - 10,
        2 * np.random.random(mna.ni) - 1,
        np.random.random(mna.nh)])    

    stuff = mna.simulate_be(10.0, 0.1, 
                            x0=x0,
                            vars=vars)

    plt.plot(stuff)
    plt.savefig('vcdcg.png')


def test_and_gate():
    solc_ = solc.Solc()
    
    a = solc_.add_variable('a')
    b = solc_.add_variable('b')
    ab = solc_.and_gate('ab', a, b)

    solc_.set_on(b)

    mna = solc_.circuit.assemble_mna_equation()

    stuff = mna.simulate_be(0.1, 0.001,
                            x0=None,
                            vars=[Var(solc_.variables[a], 'v'),
                                  Var(solc_.variables[b], 'v'),
                                  Var(solc_.variables[ab], 'v')])

    plt.plot(stuff)
    plt.savefig('and_dyn.png')


def test_factorization():
    circ = solc.Solc()

    a0 = circ.add_variable('a0')
    a1 = circ.add_variable('a1')
    b0 = circ.add_variable('b0')
    b1 = circ.add_variable('b1')

    #        A1   A0
    #        B1   B0
    # --------------
    #      A1B0 A0B0
    # A1B1 A0B1

    ab0 = circ.and_gate('ab0', a0, b0)
    ab1_1 = circ.and_gate('ab1_1', a0, b1)
    ab1_2 = circ.and_gate('ab1_2', a1, b0)
    ab1 = circ.xor_gate('ab1', ab1_1, ab1_2)
    ab1_carry = circ.and_gate('ab1_carry', ab1_1, ab1_2)
    ab2_1 = circ.and_gate('ab2_1', a1, b1)
    ab2 = circ.xor_gate('ab2', ab2_1, ab1_carry)
    ab3 = circ.and_gate('ab3', ab2_1, ab1_carry)

    # set output to 4 = 0 1 0 0
    # set output to 9 = 1 0 0 1
    # set output to 7 = 0 1 1 1
    n = 6
    
    if n & 8:
        circ.set_on(ab3)
    else:
        circ.set_off(ab3)
    if n & 4:
        circ.set_on(ab2)
    else:
        circ.set_off(ab2)
    if n & 2:
        circ.set_on(ab1)
    else:
        circ.set_off(ab1)
    if n & 1:
        circ.set_on(ab0)
    else:
        circ.set_off(ab0)


    mna = circ.circuit.assemble_mna_equation()

    stuff = mna.simulate_be(4e-4, 1e-4,
                            x0=None,
                            vars=[Var(circ.variables[a0], 'v'),
                                  Var(circ.variables[a1], 'v'),
                                  Var(circ.variables[b0], 'v'),
                                  Var(circ.variables[b1], 'v'),
                                  Var(circ.vcdcg[a0].id, 'i'),
                                  Var(circ.vcdcg[a1].id, 'i'),
                                  Var(circ.vcdcg[b0].id, 'i'),
                                  Var(circ.vcdcg[b1].id, 'i'),
                                  Var(circ.vcdcg[a0].id, 'h'),
                                  Var(circ.vcdcg[a1].id, 'h'),
                                  Var(circ.vcdcg[b0].id, 'h'),
                                  Var(circ.vcdcg[b1].id, 'h'),
                              ])

    a0_val = stuff[-1, 0] >= 0
    a1_val = stuff[-1, 1] >= 0
    b0_val = stuff[-1, 2] >= 0
    b1_val = stuff[-1, 3] >= 0

    print 'N = ', n
    print 'A = ', 2 * a1_val + a0_val
    print 'B = ', 2 * b1_val + b0_val

    plt.plot(stuff[9 * len(stuff) // 10:])
    plt.savefig('factorization.png')

    # crkt.draw()
    # plt.savefig('and.png')

if __name__ == '__main__':
    #test_vcdcgs()
    #test_and_gate()
    test_factorization()
