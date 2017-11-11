import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from anarres.circuit import Circuit, Var
from anarres.components.all import Vcdcg, Resistor
import anarres.solc as solc
import anarres.sat as sat

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
        #vars.append(Var(vcdcg.id, 'i'))
        vars.append(Var(vcdcg.id, 'h'))

    mna = crkt.assemble_mna_equation()
    v0 = 20 * np.random.randn(mna.nv)
    x0 = np.concatenate([
        v0, #np.zeros(mna.nv),
        v0[:-1], # np.zeros(mna.ni),
        np.random.random(mna.nh)])    

    # x0 = 1e-2 * np.random.randn(mna.n)
    # x0[mna.nv:mna.nv + mna.ni] = 0
    # x0[mna.nv + mna.ni:] = 0.5 + 1e-3 * np.random.randn(mna.nh)

    stuff = mna.simulate_be(2.0, 1e-3, 
                            x0=x0,
                            vars=vars)

    plt.plot(stuff)
    plt.savefig('vcdcg.png')


def test_and_gate():
    solc_ = solc.Solc()
    
    a = solc_.add_variable('a')
    b = solc_.add_variable('b')
    ab = solc_.xor_gate(a, b, name='ab')

    solc_.set_on(ab)
    solc_.set_on(a)

    mna = solc_.circuit.assemble_mna_equation()

    stuff = mna.simulate_be(0.1, 0.001,
                            x0=None,
                            vars=[Var(solc_.variables[a], 'v'),
                                  Var(solc_.variables[b], 'v'),
                                  Var(solc_.variables[ab], 'v')])

    plt.plot(stuff)
    plt.savefig('and_dyn.png')

def test_3sat():
    sat_ = sat.SoSat(3)
    sat_.add_clause(1, 2, 3)
    sat_.add_clause(1, 2, -3)
    sat_.add_clause(1, -2, 3)
    sat_.add_clause(1, -2, -3)
    sat_.add_clause(-1, 2, 3)
    sat_.add_clause(-1, 2, -3)
    sat_.add_clause(-1, -2, 3)
    #sat_.add_clause(-1, -2, -3)

    mna = sat_.circ.circuit.assemble_mna_equation()

    stuff = mna.simulate_be(0.2, 0.001,
                            x0=None,
                            vars=[Var(sat_.circ.variables['a1'], 'v'),
                                  Var(sat_.circ.variables['a2'], 'v'),
                                  Var(sat_.circ.variables['a3'], 'v')])
                                  

    plt.plot(stuff)
    plt.savefig('3sat.png')

def test_addition():
    circ = solc.Solc()

    a = circ.add_variable('a')
    b = circ.add_variable('b')

    aplusb0 = circ.xor_gate(a, b, name='aplusb0')
    aplusb1 = circ.and_gate(a, b, name='aplusb1')

    # set output to 2 = 1 0
    n = 1
    
    if n & 2:
        circ.set_on(aplusb1)
    else:
        circ.set_off(aplusb1)
    if n & 1:
        circ.set_on(aplusb0)
    else:
        circ.set_off(aplusb0)


    mna = circ.circuit.assemble_mna_equation()

    x0 = 1e0 * np.random.randn(mna.n)
    x0[mna.nv + mna.ni:] = 0.5 + 1e-3 * np.random.randn(mna.nh)

    # x0[mna.var_lut[Var(circ.variables[a], 'v')]] = 10.0
    # x0[mna.var_lut[Var(circ.variables[b], 'v')]] = -10.0

    # x0[mna.var_lut[Var(circ.vcdcg[a].id, 'i')]] = -10.0
    # x0[mna.var_lut[Var(circ.vcdcg[b].id, 'i')]] = 10.0

    # x0[mna.var_lut[Var(circ.vcdcg[a].id, 'h')]] = 0.5
    # x0[mna.var_lut[Var(circ.vcdcg[b].id, 'h')]] = 0.5


    stuff = mna.simulate_be(2e0, 1e-2,
                            x0=x0,
                            vars=[Var(circ.variables[a], 'v'),
                                  Var(circ.variables[b], 'v'),
                                  Var(circ.vcdcg[a].id, 'i'),
                                  Var(circ.vcdcg[b].id, 'i'),
                                  Var(circ.vcdcg[a].id, 'h'),
                                  Var(circ.vcdcg[b].id, 'h'),
                              ])

    a_val = int(stuff[-1, 0] >= 0)
    b_val = int(stuff[-1, 1] >= 0)

    print 'N = ', n, ' = A + B'
    print 'A = ', a_val
    print 'B = ', b_val

    plt.plot(stuff)
    plt.savefig('addition.png')


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

    ab0 = circ.and_gate(a0, b0, name='ab0')
    ab1_1 = circ.and_gate(a0, b1, name='ab1_1')
    ab1_2 = circ.and_gate(a1, b0, name='ab1_2')
    ab1 = circ.xor_gate(ab1_1, ab1_2, name='ab1')
    ab1_carry = circ.and_gate(ab1_1, ab1_2, name='ab1_carry')
    ab2_1 = circ.and_gate(a1, b1, name='ab2_1')
    ab2 = circ.xor_gate(ab2_1, ab1_carry, name='ab2')
    ab3 = circ.and_gate(ab2_1, ab1_carry, name='ab3')

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

    x0 = 1e-2 * np.random.randn(mna.n)
    x0[mna.nv + mna.ni:] = 0.5

    # x0 = np.zeros(mna.n)
    # x0[mna.var_lut[Var(circ.variables[a0], 'v')]] = 10.0
    # x0[mna.var_lut[Var(circ.variables[a1], 'v')]] = 10.0
    # x0[mna.var_lut[Var(circ.variables[b0], 'v')]] = -10.0
    # x0[mna.var_lut[Var(circ.variables[b1], 'v')]] = 10.0

    # x0[mna.var_lut[Var(circ.vcdcg[a0].id, 'i')]] = -10.0
    # x0[mna.var_lut[Var(circ.vcdcg[a1].id, 'i')]] = -10.0
    # x0[mna.var_lut[Var(circ.vcdcg[b0].id, 'i')]] = 10.0
    # x0[mna.var_lut[Var(circ.vcdcg[b1].id, 'i')]] = -10.0

    # x0[mna.var_lut[Var(circ.vcdcg[a0].id, 'h')]] = 0.5
    # x0[mna.var_lut[Var(circ.vcdcg[a1].id, 'h')]] = 0.5
    # x0[mna.var_lut[Var(circ.vcdcg[b0].id, 'h')]] = 0.5
    # x0[mna.var_lut[Var(circ.vcdcg[b1].id, 'h')]] = 0.5

    stuff = mna.simulate_be(2, 1e-1,
                            x0=x0,
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

def test_factorization3():
    import anarres.factorization as fact

    n = 2
    N = 9

    mult = fact.SoMultiplier(n)
    mult.set_answer(N)

    mna = mult.circ.circuit.assemble_mna_equation()

    x0 = 1e-2 * np.random.randn(mna.n)
    x0[mna.nv + mna.ni:] = 0.5 + 1e-3 * np.random.randn(mna.nh)

    (a0, a1), (b0, b1) = mult.get_input_vars()
    circ = mult.circ

    def hook(xs):
        a0_val, a1_val, b0_val, b1_val = [int(x >= 0) for x in xs[:4]]

        print
        print 'N = ', N
        print 'A = ', 2 * a1_val + a0_val
        print 'B = ', 2 * b1_val + b0_val

    #x0 = mna.solve_dc()

    stuff = mna.simulate_be(2, 1e-3,
                            x0=x0,
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
                              ],
                            hook=hook)

    a0_val = stuff[-1, 0] >= 0
    a1_val = stuff[-1, 1] >= 0
    b0_val = stuff[-1, 2] >= 0
    b1_val = stuff[-1, 3] >= 0

    plt.plot(stuff[9 * len(stuff) // 10:])
    plt.savefig('factorization.png')

    # crkt.draw()
    # plt.savefig('and.png')

if __name__ == '__main__':
    import cProfile

    #cProfile.run('test_vcdcgs()', 'stats')
    #test_vcdcgs()
    #test_and_gate()
    #test_addition()
    #test_factorization3()
    
    test_3sat()
