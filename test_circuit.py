import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import circuit


def test_circuit_resistance():
    crkt = circuit.Circuit()
    battery = circuit.VoltageSource(crkt, 9)  # 9 volts
    resistor = circuit.Resistor(crkt, 10)  # 1 ohm    
    battery.negative = resistor.positive
    battery.positive = resistor.negative = crkt._ground
    print crkt.graph.edges()
    mna = crkt.assemble_mna_equation()
    stuff = mna.simulate(10.0, 0.1)

    crkt.draw()
    plt.savefig('crkt1.png')

    
def test_circuit_capacitance():
    crkt = circuit.Circuit()
    battery = circuit.VoltageSource(crkt, 9)
    resistor = circuit.Resistor(crkt, 1)
    capacitor = circuit.Capacitor(crkt, 1e-1)
    battery.positive = resistor.positive
    resistor.negative = capacitor.positive
    capacitor.negative = battery.negative = crkt._ground

    mna = crkt.assemble_mna_equation()
    stuff = mna.simulate_be(1e0, 1e-3, x0=[0, 9, 0], vars=['node_7', 'node_3'])

    plt.plot(stuff)
    plt.savefig('cap.png')

    plt.figure()
    crkt.draw()
    plt.savefig('crkt2.png')

def test_circuit_parallel_resistor_capacitor():
    crkt = circuit.Circuit()
    battery = circuit.VoltageSource(crkt, 9)
    resistor = circuit.Resistor(crkt, 1)
    capacitor = circuit.Capacitor(crkt, 1e-1)  # 1 microFarad
    battery.negative = resistor.positive = capacitor.positive
    battery.positive = resistor.negative = capacitor.negative = crkt._ground

    mna = crkt.assemble_mna_equation()
    stuff = mna.simulate_be(1e0, 1e-3, x0=[9, 0], vars=['V_2'])

    print mna.rev_idict

    plt.plot(stuff)
    plt.savefig('cap_par.png')

    plt.figure()
    crkt.draw()
    plt.savefig('crkt3.png')


if __name__ == '__main__':
    #test_circuit_resistance()
    #test_circuit_capacitance()
    test_circuit_parallel_resistor_capacitor()
