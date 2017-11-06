import circuit2

class Solc(object):
    def __init__(self):
        self.circuit = circuit2.Circuit()
        self.variables = dict()
        self.vcdcg = dict()

    def add_variable(self, name):
        aux = circuit2.Vcdcg(self.circuit)
        aux.negative = self.circuit._ground
        self.variables[name] = aux.positive
        self.vcdcg[name] = aux
        shunt = circuit2.Resistor(self.circuit, 1e9)
        shunt.negative = self.circuit._ground
        shunt.positive = self.variables[name]
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
        hardwired_voltage.negative = self.circuit._ground
        hardwired_voltage.positive = self.variables[name]

    def set_off(self, name):
        hardwired_voltage = circuit2.VoltageSource(self.circuit, -10)
        hardwired_voltage.negative = self.circuit._ground
        hardwired_voltage.positive = self.variables[name]
