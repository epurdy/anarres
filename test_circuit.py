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
    print stuff

    
def test_circuit_capacitance():
    crkt = circuit.Circuit()
    battery = circuit.VoltageSource(crkt, 9)
    resistor = circuit.Resistor(crkt, 1)
    capacitor = circuit.Capacitor(crkt, 1e-6)  # 1 microFarad
    battery.negative = resistor.positive
    resistor.negative = capacitor.positive
    capacitor.negative = battery.positive = crkt._ground

    mna = crkt.assemble_mna_equation()
    stuff = mna.simulate(0.01, 0.001)
    print mna.rev_vdict
    print mna.rev_idict

def test_circuit_parallel_resistor_capacitor():
    crkt = circuit.Circuit()
    battery = crky.VoltageSource(9)
    resistor = crkt.Resistor(1)
    capacitor = crky.Capacitor(1e-6)  # 1 microFarad
    battery.negative = resistor.positive = capacitor.positive
    battery.positive = resistor.negative = capacitor.negative


if __name__ == '__main__':
    #test_circuit_resistance()
    test_circuit_capacitance()
