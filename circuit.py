from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


NEXT_ID = 0
def get_next_id(typ):
    global NEXT_ID
    rv = '{}_{}'.format(typ, NEXT_ID)
    NEXT_ID += 1
    print('ID', rv)
    return rv


class MnaEquation(object):
    def __init__(self, A, B, c, rev_vdict, rev_idict):
        self.A = A
        self.B = B
        self.c = c
        self.rev_vdict = rev_vdict
        self.rev_idict = rev_idict
        self.nv = len(self.rev_vdict)
        self.ni = len(self.rev_idict)

    def simulate(self, tmax, dt):
        x, _, _, _ = np.linalg.lstsq(self.A, self.c)
        xs = [x]
        t = 0.0
        while t < tmax:            
            dx, _, _, _ = np.linalg.lstsq(self.B, self.c - self.A.dot(x))
            print dx
            x += dt * dx
            xs.append(x)
            t += dt

        return np.array(xs)

class Circuit(object):
    def __init__(self):
        self.id = get_next_id('crkt')
        self.graph = nx.MultiDiGraph()
        self.components = []
        self.item_lut = dict()
        self.node_lut = dict()
        
        # add ground
        self._ground = get_next_id('node')
        self.add_node(self._ground)
        
    def add_node(self, node):
        self.graph.add_node(node)
        self.node_lut[node] = node
        
    def add_component(self, comp, nodea=None, nodeb=None):
        if nodea is not None:
            self.add_node(nodea)

        if nodeb is not None:
            self.add_node(nodeb)

        if nodea is not None and nodeb is not None:
            self.graph.add_edge(nodea, nodeb, data=comp.id)

        self.components.append(comp)
        self.item_lut[comp.id] = comp
        
    def merge(self, nodea, nodeb):
        nodea = self.node_lut[nodea]
        nodeb = self.node_lut[nodeb]

        if nodea > nodeb:
            nodea, nodeb = nodeb, nodea

        self.node_lut[nodeb] = nodea

        self.graph = nx.contracted_nodes(self.graph, nodea, nodeb, self_loops=False)

    def assemble_mna_equation(self):
        vvars = list(x for x in self.graph.nodes() if x != self._ground)
        ivars = list(x.id for x in self.components if x.has_current_variable)
        nv = len(vvars)
        ni = len(ivars)
        rev_vdict = {id:num for num, id in enumerate(vvars)}
        rev_idict = {id:num + nv for num, id in enumerate(ivars)}

        # A x + B x' = c        
        A = np.zeros((nv + ni, nv + ni))
        B = np.zeros((nv + ni, nv + ni))
        c = np.zeros(nv + ni)

        def resolve_id_var(id, var):
            if var == 'v':
                id = self.node_lut[id]
                if id == self._ground:
                    return None
                else:
                    return rev_vdict[id]
            elif var == 'i':
                return rev_idict[id]
            else:
                assert False, 'Variable must be either v or i'

        def use_stamp1(mat, iid, ivar, val):
            i = resolve_id_var(iid, ivar)
            if i is not None:
                mat[i] += val

        def use_stamp2(mat, iid, ivar, jid, jvar, val):
            i = resolve_id_var(iid, ivar)
            j = resolve_id_var(jid, jvar)
            if i is not None and j is not None:
                mat[i, j] += val

        for component in self.components:
            for stuff in component.Astamp():
                use_stamp2(A, *stuff)

            for stuff in component.Bstamp():
                use_stamp2(B, *stuff)

            for stuff in component.cstamp():
                use_stamp1(c, *stuff)

        return MnaEquation(A=A, B=B, c=c, 
                           rev_vdict=rev_vdict, 
                           rev_idict=rev_idict)

        
class CircuitComponent(object):
    __metaclass__ = ABCMeta

    def __init__(self, crkt):
        self._circuit = crkt

    def Astamp(self):
        return []

    def Bstamp(self):
        return []

    def cstamp(self):
        return []

    @property
    @abstractmethod
    def has_current_variable(self):
        pass

class TwoNodeCircuitComponent(CircuitComponent):
    def __init__(self, crkt):
        super(TwoNodeCircuitComponent, self).__init__(crkt)
        self._positive = get_next_id('node')
        self._negative = get_next_id('node')
        crkt.add_component(self, self._positive, self._negative)
        
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

    
class Resistor(TwoNodeCircuitComponent):
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor, self).__init__(crkt)
        self.resistance = float(resistance)

    @property
    def has_current_variable(self):
        return False

    def Astamp(self):
        return [
            (self.positive, 'v', self.positive, 'v',  1.0 / self.resistance),
            (self.positive, 'v', self.negative, 'v', -1.0 / self.resistance),
            (self.negative, 'v', self.positive, 'v', -1.0 / self.resistance),
            (self.negative, 'v', self.negative, 'v',  1.0 / self.resistance),
            ]


class Capacitor(TwoNodeCircuitComponent):
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor, self).__init__(crkt)
        self.capacitance = float(capacitance)

    @property
    def has_current_variable(self):
        return False

    def Bstamp(self):
        return [
            (self.positive, 'v', self.positive, 'v',  self.capacitance),
            (self.positive, 'v', self.negative, 'v', -self.capacitance),
            (self.negative, 'v', self.positive, 'v', -self.capacitance),
            (self.negative, 'v', self.negative, 'v',  self.capacitance),
            ]

    
class VoltageSource(TwoNodeCircuitComponent):
    def __init__(self, crkt, voltage):
        self.id = get_next_id('V')
        super(VoltageSource, self).__init__(crkt)
        self.voltage = voltage

    @property
    def has_current_variable(self):
        return True

    def Astamp(self):
        return [
            (self.positive, 'v', self.id, 'i',   1.0),
            (self.negative, 'v', self.id, 'i',  -1.0),
            (self.id, 'i', self.positive, 'v',   1.0),
            (self.id, 'i', self.negative, 'v',  -1.0),
        ]

    def cstamp(self):
        return [
            (self.id, 'i', self.voltage),
        ]


# class Memristor(TwoNodeCircuitComponent):
#     def M(self, x):
#         return self.R_on * (1 - x) + self.R_off * x

        
# class Vcvg(CircuitComponent):
#     pass


# class Vcdcg(CircuitComponent):
#     pass

    
# class DynamicCorrectionModule(CircuitComponent):
#     pass


# class SoGate(CircuitComponent):
#     """
#     v1--+   +--v2
#         |   |
#  G--DCM-+   +-DCM--G
#         |   |
#         R   R
#         +_+_+
#           |
#   G--DCM--+
#           |
#           vo
#     """
#     def __init__(self, v1, v2):
#         self.v1 = v1
#         self.v2 = v2
        
# class SoAnd(SoGate):
#     pass
# class SoOr(SoGate):
#     pass
# class SoXor(SoGate):
#     pass
        

    
