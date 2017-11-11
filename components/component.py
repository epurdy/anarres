from abc import ABCMeta, abstractmethod

from util import get_next_id
        
class CircuitComponent(object):
    __metaclass__ = ABCMeta

    def __init__(self, crkt):
        self._circuit = crkt

    def Astamp(self, static=False):
        return []

    def Bstamp(self):
        return []

    def cstamp(self, static=False):
        return []

    @abstractmethod
    def has_current_variable(self):
        pass

    @abstractmethod
    def has_hidden_variable(self):
        pass


class TwoNodeCircuitComponent(CircuitComponent):
    def __init__(self, crkt):
        super(TwoNodeCircuitComponent, self).__init__(crkt)
        self._positive = get_next_id('node')
        self._negative = get_next_id('node')
        crkt.add_component(self, self._positive, self._negative)

    def normalize(self):
        redo = False
        pos = self._circuit.node_lut[self._positive]
        if self._positive != pos:
            self._positive = pos
            redo = True

        neg = self._circuit.node_lut[self._negative]
        if self._negative != neg:
            self._negative = neg
            redo = True

        return redo
        
    @property
    def positive(self):
        return self._positive

    @positive.setter
    def positive(self, other):
        self._circuit.merge(self._positive, other)

    @property
    def negative(self):
        return self._negative

    @negative.setter
    def negative(self, other):
        self._circuit.merge(self._negative, other)

class ThreeNodeCircuitComponent(CircuitComponent):
    def __init__(self, crkt):
        super(ThreeNodeCircuitComponent, self).__init__(crkt)
        self._node1 = get_next_id('node')
        self._node2 = get_next_id('node')
        self._node3 = get_next_id('node')
        # note that we do not add the component to the circuit's edge set!
        crkt.add_node(self._node1)
        crkt.add_node(self._node2)
        crkt.add_node(self._node3)
        
    @property
    def node1(self):
        return self._node1

    @node1.setter
    def node1(self, other):
        self._circuit.merge(self._node1, other)

    @property
    def node2(self):
        return self._node2

    @node2.setter
    def node2(self, other):
        self._circuit.merge(self._node2, other)

    @property
    def node3(self):
        return self._node3

    @node3.setter
    def node3(self, other):
        self._circuit.merge(self._node3, other)

    def normalize(self):
        redo = False

        node1 = self._circuit.node_lut[self._node1]
        if self._node1 != node1:
            self._node1 = node1
            redo = True

        node2 = self._circuit.node_lut[self._node2]
        if self._node2 != node2:
            self._node2 = node2
            redo = True

        node3 = self._circuit.node_lut[self._node3]
        if self._node3 != node3:
            self._node3 = node3
            redo = True

        return redo
