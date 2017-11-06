import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from circuit import Circuit, Var
from components.all import Resistor, Capacitor, VoltageSource

def test_circuit_resistance():
    circuit = Circuit()
    battery = VoltageSource(circuit, 9)  # 9 volts
    resistor = Resistor(circuit, 10)  # 1 ohm    
    battery.negative = resistor.positive
    battery.positive = resistor.negative = circuit._ground
    print circuit.graph.edges()
    mna = circuit.assemble_mna_equation()
    stuff = mna.simulate(10.0, 0.1)

    circuit.draw()
    plt.savefig('crkt1.png')

    
def test_circuit_capacitance():
    circuit = Circuit()
    battery = VoltageSource(circuit, 9)
    resistor = Resistor(circuit, 1)
    capacitor = Capacitor(circuit, 1e-1)
    battery.positive = resistor.positive
    resistor.negative = capacitor.positive
    capacitor.negative = battery.negative = circuit._ground

    mna = circuit.assemble_mna_equation()
    stuff = mna.simulate_be(1e0, 1e-3, 
                            vars=[Var(resistor.negative, 'v'),
                                  Var(resistor.positive, 'v')])

    plt.plot(stuff)
    plt.savefig('cap.png')

    plt.figure()
    circuit.draw()
    plt.savefig('crkt2.png')

def test_circuit_parallel_resistor_capacitor():
    circuit = Circuit()
    battery = VoltageSource(circuit, 9)
    resistor = Resistor(circuit, 1)
    capacitor = Capacitor(circuit, 1e-1)  # 1 microFarad
    battery.negative = resistor.positive = capacitor.positive
    battery.positive = resistor.negative = capacitor.negative = circuit._ground

    mna = circuit.assemble_mna_equation()
    stuff = mna.simulate_be(1e0, 1e-3, 
                            vars=[Var(battery.id, 'i')])

    print mna.rev_idict

    plt.plot(stuff)
    plt.savefig('cap_par.png')

    plt.figure()
    circuit.draw()
    plt.savefig('crkt3.png')


if __name__ == '__main__':
    #test_circuit_resistance()
    #test_circuit_capacitance()
    test_circuit_parallel_resistor_capacitor()
