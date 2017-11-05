from abc import ABCMeta, abstractmethod
from collections import defaultdict
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

def indent_str(s, n=2):
    lines = s.split('\n')
    lines = [' ' * n + x for x in lines]
    return '\n'.join(lines)

class MnaEquation(object):
    def __init__(self, A, B, c, rev_vdict, rev_idict):
        self.A = A
        self.B = B
        self.c = c
        self.rev_vdict = rev_vdict
        self.rev_idict = rev_idict
        self.nv = len(self.rev_vdict)
        self.ni = len(self.rev_idict)

        self.rev_dict = {}
        for k in self.rev_vdict:
            self.rev_dict[k] = self.rev_vdict[k]
        for k in self.rev_idict:
            self.rev_dict[k] = self.rev_idict[k]

    def __repr__(self):
        Astr = indent_str(str(self.A), 14)
        Bstr = indent_str(str(self.B), 14)
        Binvstr = indent_str(str(np.linalg.pinv(self.B)), 14)
        cstr = indent_str(str(self.c), 14)
        return """
        A x + B x' = c
          A = 
{}
          B = 
{}
          Binv = 
{}
          c = 
{}
        """.format(Astr, Bstr, Binvstr, cstr)

    def simulate_be(self, tmax, dt, x0=None, vars=None):
        if x0 is None:
            x, res, rank, _ = np.linalg.lstsq(self.A, self.c)
            print res, rank
            print x
            print 'residual', self.c - self.A.dot(x)
            print 'residual', self.c - self.A.dot([0, 9, -9, 9, 9])
        else:
            x = np.array(x0)

        if vars is None:
            indices = range(self.nv + self.ni)
        else:
            indices = [self.rev_dict[k] for k in vars]

        xs = [x[indices]]
        t = 0.0
        n = self.nv + self.ni
        mat = self.A + self.B / dt

        while t < tmax:
            x, _, _, _ = np.linalg.lstsq(mat, self.c + (1.0/dt) * self.B.dot(x))
            xs.append(list(x[indices]))
            t += dt
            #print x
            #print xdx

        xs = np.array(xs)
        print xs.shape
        print xs

        return xs

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

    def draw(self):
        pos = nx.drawing.layout.spring_layout(self.graph)
        edge_labels = defaultdict(lambda: '')
        for edge in self.graph.edges():
            edge_data = self.graph.get_edge_data(*edge)
            for idx in edge_data:
                edge_labels[edge] += edge_data[idx]['data']
        nx.draw_networkx(self.graph, pos, with_labels=True)
        nx.draw_networkx_edge_labels(self.graph, pos, edge_labels=edge_labels)
        
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
    """A resistor without an explicit current variable.
    """
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


class Resistor2(TwoNodeCircuitComponent):
    """A resistor with an explicit current variable.
    """
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor2, self).__init__(crkt)
        self.resistance = float(resistance)

    @property
    def has_current_variable(self):
        return True

    def Astamp(self):
        return [
            (self.positive, 'v', self.id, 'i',  1),
            (self.negative, 'v', self.id, 'i', -1),
            (self.id, 'i', self.positive, 'v',  1),
            (self.id, 'i', self.negative, 'v', -1),
            (self.id, 'i', self.id, 'i',  -self.resistance),
        ]


class Capacitor(TwoNodeCircuitComponent):
    """A capacitor without an explicit current variable.
    """
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


class Capacitor2(TwoNodeCircuitComponent):
    """A capacitor with an explicit current variable.
    """
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor2, self).__init__(crkt)
        self.capacitance = float(capacitance)

    @property
    def has_current_variable(self):
        return True

    def Astamp(self):
        return [
            # these do not appear in Najm but seem to be necessary
            (self.positive, 'v', self.id, 'i',  1),
            (self.negative, 'v', self.id, 'i', -1),
            (self.id, 'i', self.positive, 'v',  1),
            (self.id, 'i', self.negative, 'v', -1),

            # this one does appear in Najm
            (self.id, 'i', self.id, 'i',  1),
        ]

    def Bstamp(self):
        return [
            (self.id, 'i', self.positive, 'v',  self.capacitance),
            (self.id, 'i', self.negative, 'v', -self.capacitance),
        ]

    
class VoltageSource(TwoNodeCircuitComponent):
    """Voltage sources are required to have an explicit current variable."""
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
        

    
