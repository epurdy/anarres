"""
A higher-level API, with which one can build self-organizing versions of standard logic circuits.
"""

from anarres.circuit import Circuit
from anarres.components.all import Resistor, VoltageSource, Vcdcg, SoAnd, SoOr, SoXor

class Solc(object):
    def __init__(self):
        self.circuit = Circuit()
        self.variables = dict()
        self.vcdcg = dict()

    def add_variable(self, name):
        aux = Vcdcg(self.circuit)
        aux.negative = self.circuit._ground
        self.variables[name] = aux.positive
        self.vcdcg[name] = aux
        shunt = Resistor(self.circuit, 1e9)
        shunt.negative = self.circuit._ground
        shunt.positive = self.variables[name]
        return name

    def and_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = SoAnd(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def or_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = SoOr(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def xor_gate(self, name, in1, in2):
        name = self.add_variable(name)
        gate = SoXor(self.circuit)
        gate.node1 = self.variables[in1]
        gate.node2 = self.variables[in2]
        gate.node3 = self.variables[name]
        return name

    def set_on(self, name):
        hardwired_voltage = VoltageSource(self.circuit, 10)        
        hardwired_voltage.negative = self.circuit._ground
        hardwired_voltage.positive = self.variables[name]

    def set_off(self, name):
        hardwired_voltage = VoltageSource(self.circuit, -10)
        hardwired_voltage.negative = self.circuit._ground
        hardwired_voltage.positive = self.variables[name]
