from abc import ABCMeta, abstractmethod
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import scipy

import dcm
import vcdcg

"""circuit.py can handle linear systems with no time dependence and no extra variables. 

this should handle non-linear systems with extra variables. I don't think we
need explicit time dependence.

the MNA equation goes from being Ax + Bx' = c with A, B, c constant numpy
arrays to A,B,c being collections of stamps. A stamp is an i, j pair and a
function to call to get the value of the thing.

Need some way to registre variables so that the MNA equation knows what its
variables are.

Also, it seems like the symmetry of the is and the js breaks down... i is
really indexing into equations and j is indexing into variables; but, the
variables happen to be in 1-1 corr with the equations for now... not sure why
or if this is a general phenomenon.

We need some way of grouping variables or sth so that we can have numpy handle
things efficiently.

"""

R_ON = 0.1
R_OFF = 0.2


NEXT_ID = 0
def get_next_id(typ):
    global NEXT_ID
    rv = '%s_%06d' % (typ, NEXT_ID)
    NEXT_ID += 1
    print('ID', rv)
    return rv

class Var(object):
    def __init__(self, id, quantity):
        self.id = id
        self.quantity = quantity

    def __repr__(self):
        return '{}.{}'.format(self.id, self.quantity)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return str(self) != str(other)

    def __hash__(self):
        return hash(str(self))

class FuncStamp1(object):
    def __init__(self, ivar, fn, input_vars):
        self.ivar = ivar
        self.fn = fn
        self.input_vars = input_vars

    def __repr__(self):
        return 'FuncStamp1({}, {})'.format(self.ivar, self.fn.__name__)

    def __call__(self, x, var_lut):
        inputs = [x[var_lut[var]] 
                  if var_lut[var] is not None
                  else 0.0
                  for var in self.input_vars
              ]
        return self.fn(*inputs)

class FuncStamp2(FuncStamp1):
    def __init__(self, ivar, jvar, fn, input_vars):
        self.ivar = ivar
        self.jvar = jvar
        self.fn = fn
        self.input_vars = input_vars

    def __repr__(self):
        return 'FuncStamp2({}, {}, {})'.format(self.ivar, self.jvar, self.fn.__name__)

class ConstStamp1(object):
    def __init__(self, ivar, val):
        self.ivar = ivar
        self.val = val

    def __repr__(self):
        return 'ConstStamp1({}, {})'.format(self.ivar, self.val)

    def __call__(self, x, var_lut):
        return self.val

class ConstStamp2(ConstStamp1):
    def __init__(self, ivar, jvar, val):
        self.ivar = ivar
        self.jvar = jvar
        self.val = val

    def __repr__(self):
        return 'ConstStamp2({}, {}, {})'.format(self.ivar, self.jvar, self.val)


class MnaEquation(object):
    def __init__(self, A, B, c, rev_vdict, rev_idict, rev_hdict):
        self.A = A
        self.B = B
        self.c = c
        self.rev_vdict = rev_vdict
        self.rev_idict = rev_idict
        self.rev_hdict = rev_hdict
        self.nv = len(self.rev_vdict)
        self.ni = len(self.rev_idict)
        self.nh = len(self.rev_hdict)
        self.n = self.nv + self.ni + self.nh

        self.var_lut = {}
        for k in self.rev_vdict:
            self.var_lut[Var(k, 'v')] = self.rev_vdict[k]
        for k in self.rev_idict:
            self.var_lut[Var(k, 'i')] = self.rev_idict[k]
        for k in self.rev_hdict:
            self.var_lut[Var(k, 'h')] = self.rev_hdict[k]

    def mat(self, stamps, x):
        #rv = np.zeros((self.n, self.n))
        rv = scipy.sparse.dok_matrix((self.n, self.n), dtype=np.float32)
        for stamp in stamps:
            i = self.var_lut[stamp.ivar]
            j = self.var_lut[stamp.jvar]

            # don't include ground
            if i is None or j is None:
                continue

            val = stamp(x, self.var_lut)
            rv[i, j] += val
        return rv.tocoo()

    def Amat(self, x):
        return self.mat(self.A, x)

    def Bmat(self, x):
        return self.mat(self.B, x)

    def cmat(self, x):
        rv = np.zeros(self.n)
        #rv = scipy.sparse.dok_matrix(self.n, dtype=np.float32)
        for stamp in self.c:
            i = self.var_lut[stamp.ivar]

            # don't include ground
            if i is None:
                continue

            val = stamp(x, self.var_lut)
            rv[i] += val
        return rv

    def simulate_be(self, tmax, dt, x0=None, vars=None):
        if x0 is None:
            xm1 = np.zeros(self.n)
            x, res, rank, _ = np.linalg.lstsq(self.Amat(xm1), self.cmat(xm1))
        else:
            x = np.array(x0)

        if vars is None:
            indices = range(self.n)
        else:
            indices = [self.var_lut[k] for k in vars]

        xs = [x[indices]]
        t = 0.0

        thresh = 10.0
        while t < tmax:
            print t

            oldx = x
            def f(x):
                # A(x) x + B(x) x' = c(x)
                # A(x) x + B(x) * (x - oldx)/dt = c(x)
                # (A(x) + B(x)/dt) x = c(x) + B(x)/dt oldx
                A = self.Amat(x)
                B = self.Bmat(x)
                c = self.cmat(x)
                #return (A + B/dt).dot(x) - c - (B/dt).dot(oldx)
                return A.dot(x) + B.dot(x - oldx) / dt - c

            x, info, ier, mesg = scipy.optimize.fsolve(f, oldx, full_output=True)
#                                                       maxfev=3 * (self.n + 1))
            print np.linalg.norm(f(x))
            print mesg

            #x += 1e-6 * np.random.randn(self.n)
            
            xs.append(list(x[indices]))
            print zip(vars, xs[-1])
            t += dt

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
        graph = self.graph.copy()
        graph.remove_node(self._ground)
        # pos = nx.drawing.layout.spring_layout(self.graph)
        # pos = nx.drawing.layout.circular_layout(self.graph)
        pos = nx.drawing.layout.shell_layout(self.graph)
        # pos = nx.drawing.layout.spectral_layout(graph)
        edge_labels = defaultdict(lambda: '')
        for edge in graph.edges():
            edge_data = graph.get_edge_data(*edge)
            for idx in edge_data:
                edge_labels[edge] += edge_data[idx]['data']
        nx.draw_networkx(graph, pos, with_labels=True)
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
        
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

    def normalize(self):
        redo = False
        for comp in self.components:
            redo = comp.normalize() or redo
        if redo:
            self.normalize()
        
    def merge(self, nodea, nodeb):
        nodea = self.node_lut[nodea]
        nodeb = self.node_lut[nodeb]

        if nodea > nodeb:
            nodea, nodeb = nodeb, nodea

        self.node_lut[nodeb] = nodea

        self.graph = nx.contracted_nodes(self.graph, nodea, nodeb, self_loops=False)

    def assemble_mna_equation(self):
        self.normalize()

        vvars = [x for x in self.graph.nodes() if x != self._ground]
        ivars = [x.id for x in self.components if x.has_current_variable()]
        hvars = [x.id for x in self.components if x.has_hidden_variable()]
        nv = len(vvars)
        ni = len(ivars)
        nh = len(hvars)
        rev_vdict = {id:num for num, id in enumerate(vvars)}
        rev_idict = {id:num + nv for num, id in enumerate(ivars)}
        rev_hdict = {id:num + nv + ni for num, id in enumerate(hvars)}

        rev_vdict[self._ground] = None

        # A x + B x' = c        
        A = []
        B = []
        c = []

        for component in self.components:
            A.extend(component.Astamp())
            B.extend(component.Bstamp())
            c.extend(component.cstamp())

        return MnaEquation(A=A, B=B, c=c, 
                           rev_vdict=rev_vdict, 
                           rev_idict=rev_idict,
                           rev_hdict=rev_hdict)

        
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

    
class Resistor(TwoNodeCircuitComponent):
    """A resistor without an explicit current variable.
    """
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor, self).__init__(crkt)
        self.resistance = float(resistance)

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        return [
            ConstStamp2(posvar, posvar,  1.0 / self.resistance),
            ConstStamp2(posvar, negvar, -1.0 / self.resistance),
            ConstStamp2(negvar, posvar, -1.0 / self.resistance),
            ConstStamp2(negvar, negvar,  1.0 / self.resistance),
        ]


class Resistor2(TwoNodeCircuitComponent):
    """A resistor with an explicit current variable.
    """
    def __init__(self, crkt, resistance):
        self.id = get_next_id('R')
        super(Resistor2, self).__init__(crkt)
        self.resistance = float(resistance)

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),
            ConstStamp2(ivar, ivar, -self.resistance)
        ]


class Capacitor(TwoNodeCircuitComponent):
    """A capacitor without an explicit current variable.
    """
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor, self).__init__(crkt)
        self.capacitance = float(capacitance)

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

    def Bstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        return [
            ConstStamp2(posvar, posvar,  self.capacitance),
            ConstStamp2(posvar, negvar, -self.capacitance),
            ConstStamp2(negvar, posvar, -self.capacitance),
            ConstStamp2(negvar, negvar,  self.capacitance),
        ]


class Capacitor2(TwoNodeCircuitComponent):
    """A capacitor with an explicit current variable.
    """
    def __init__(self, crkt, capacitance):
        self.id = get_next_id('C')
        super(Capacitor2, self).__init__(crkt)
        self.capacitance = float(capacitance)

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            # these do not appear in Najm but seem to be necessary
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),

            # this one does appear in Najm
            ConstStamp2(ivar, ivar, 1),
        ]

    def Bstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(ivar, posvar,  self.capacitance),
            ConstStamp2(ivar, negvar, -self.capacitance),
        ]

    
class VoltageSource(TwoNodeCircuitComponent):
    """Voltage sources are required to have an explicit current variable."""
    def __init__(self, crkt, voltage):
        self.id = get_next_id('V')
        super(VoltageSource, self).__init__(crkt)
        self.voltage = voltage

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),
        ]

    def cstamp(self):
        ivar = Var(self.id, 'i')
        return [
            ConstStamp1(ivar, self.voltage),
        ]


class Memristor(TwoNodeCircuitComponent):
    def __init__(self, crkt, R_on, R_off, C=1e-3):
        self.id = get_next_id('M')
        super(Memristor, self).__init__(crkt)
        self.R_on = R_on
        self.R_off = R_off

        # a capacitor in parallel to model parasitic capcitance
        self.parasitic_capacitor = Capacitor(crkt, C)
        self.parasitic_capacitor.positive = self.positive
        self.parasitic_capacitor.negative = self.negative

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return True

    def M(self, x):
        return self.R_on * (1 - x) + self.R_off * x

    def Minv(self, x):
        return 1.0 / self.M(x)

    def negMinv(self, x):
        return -1.0 / self.M(x)

    def fM(self, pos, neg, h):
        i = (pos - neg) * self.Minv(h)
        return dcm.derivative_of_memristor_variables(pos - neg, h)

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp2(posvar, posvar, self.Minv,    [hvar]),
            FuncStamp2(posvar, negvar, self.negMinv, [hvar]),
            FuncStamp2(negvar, posvar, self.negMinv, [hvar]),
            FuncStamp2(negvar, negvar, self.Minv,    [hvar]),
        ]

    def Bstamp(self):
        hvar = Var(self.id, 'h')
        return [
            ConstStamp2(hvar, hvar, 1)
        ]

    def cstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp1(hvar, self.fM, [posvar, negvar, hvar])
        ]
        
class Vcvg(TwoNodeCircuitComponent):
    def __init__(self, crkt, params, nodes):
        self.id = get_next_id('VCVG')
        super(Vcvg, self).__init__(crkt)

        self.params = params
        self.nodes = nodes

    def normalize(self): 
        redo = super(Vcvg, self).normalize()       
        norm_nodes = [self._circuit.node_lut[node] for node in self.nodes]
        if norm_nodes != self.nodes:
            redo = True
            self.nodes = norm_nodes
        return redo

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return False

    def voltage_fn(self, *x):
        return self.params.dot(x + (1,))

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
            ConstStamp2(ivar, posvar,  1),
            ConstStamp2(ivar, negvar, -1),
        ]

    def cstamp(self):
        ivar = Var(self.id, 'i')
        input_vars = [Var(node, 'v') for node in self.nodes]
        return [
            FuncStamp1(ivar, self.voltage_fn, input_vars),
        ]
        

class Vcdcg(TwoNodeCircuitComponent):
    _all_vcdcg_currents = []

    @classmethod
    def all_vcdcg_currents(cls):
        return cls._all_vcdcg_currents

    def __init__(self, crkt):
        self.id = get_next_id('VCDCG')
        super(Vcdcg, self).__init__(crkt)
        self.__class__._all_vcdcg_currents.append(Var(self.id, 'i'))

    def has_current_variable(self):
        return True

    def has_hidden_variable(self):
        return True

    def deriv_i(self, posv, negv, i, h):
        v = negv - posv
        rv = vcdcg.derivative_of_dcg_currents(v, i, h, delta_s=1e-3, gamma=1.0, v_c=10.0)
        return rv

    def deriv_h(self, h, *ivals):
        # note that the ivals are out of order but that it does not matter
        ivals = np.array(ivals)
        rv = vcdcg.derivative_of_dcg_variables(ivals, h, i_min=9.0, i_max=11.0, k_s=1.0, k_i=3.0, delta_i=1e-3)
        return rv

    def Astamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        return [
            ConstStamp2(posvar, ivar,  1),
            ConstStamp2(negvar, ivar, -1),
        ]

    def Bstamp(self):
        ivar = Var(self.id, 'i')
        hvar = Var(self.id, 'h')
        return [
            ConstStamp2(ivar, ivar, 1),
            ConstStamp2(hvar, hvar, 1),
        ]

    def cstamp(self):
        posvar = Var(self.positive, 'v')
        negvar = Var(self.negative, 'v')
        ivar = Var(self.id, 'i')
        hvar = Var(self.id, 'h')
        return [
            FuncStamp1(ivar, self.deriv_i, [posvar, negvar, ivar, hvar]),
            FuncStamp1(hvar, self.deriv_h, [hvar] + self.all_vcdcg_currents())
        ]

    
class DynamicCorrectionModule(TwoNodeCircuitComponent):
    def __init__(self, crkt, params, nodes, r_on=R_ON, r_off=R_OFF):
        self.id = get_next_id('DCM')
        super(DynamicCorrectionModule, self).__init__(crkt)

        gnd = crkt._ground

        self.negative = gnd

        self.vcvg_b = Vcvg(crkt, params[0], nodes)
        self.vcvg_c = Vcvg(crkt, params[1], nodes)
        self.vcvg_d = Vcvg(crkt, params[4], nodes)
        self.vcvg_e = Vcvg(crkt, params[2], nodes)
        self.vcvg_f = Vcvg(crkt, params[3], nodes)

        self.mem_b = Memristor(crkt, r_on, r_off)
        self.mem_c = Memristor(crkt, r_on, r_off)
        self.mem_e = Memristor(crkt, r_on, r_off)
        self.mem_f = Memristor(crkt, r_on, r_off)

        self.res_d = Resistor(crkt, r_off)

        # hook negative of vcvg's to ground
        self.vcvg_b.negative = gnd
        self.vcvg_c.negative = gnd
        self.vcvg_d.negative = gnd
        self.vcvg_e.negative = gnd
        self.vcvg_f.negative = gnd

        # hook positive of vcvg's to negative of corresponding
        # memristor/resistor
        self.vcvg_b.positive = self.mem_b.negative
        self.vcvg_c.positive = self.mem_c.negative
        self.vcvg_d.positive = self.res_d.negative
        self.vcvg_e.positive = self.mem_e.negative
        self.vcvg_f.positive = self.mem_f.negative
        
        # hook positive of memristor/resistors to positive output
        self.mem_b.positive = self.positive
        self.mem_c.positive = self.positive
        self.res_d.positive = self.positive
        self.mem_e.positive = self.positive
        self.mem_f.positive = self.positive

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False


vc = 10
GATE_PARAMS = dict(
    AND=np.array(
        [[[ 0, -1,  1,      vc],  # B
          [ 1,  0,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 1,  0,  0,       0],  # F
          [ 4,  1, -3,     -vc]], # D

         [[-1,  0,  1,      vc],  # B
          [ 0,  1,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 0,  1,  0,       0],  # F
          [ 1,  4, -3,     -vc]], # D

         [[ 1,  0,  0,       0],  # B
          [ 0,  1,  0,       0],  # C
          [ 0,  0,  1,       0],  # E
          [ 2,  2, -1, -2 * vc],  # F
          [-4, -4,  7,  2 * vc]], # D
     ]),

    OR=np.array(
        [[[ 0,  0,  1,       0],  # B
          [ 1,  0,  0,       0],  # C
          [ 0, -1,  1,     -vc],  # E
          [ 1,  0,  0,       0],  # F
          [ 4,  1, -3,      vc]], # D
         
         [[ 0,  0,  1,       0],  # B
          [ 0,  1,  0,       0],  # C
          [-1,  0,  1,     -vc],  # E
          [ 0,  1,  0,       0],  # F
          [ 1,  4, -3,      vc]], # D

         [[ 0,  0,  1,       0],  # B
          [ 2,  2, -1,  2 * vc],  # C
          [ 1,  0,  0,       0],  # E
          [ 0,  1,  0,       0],  # F
          [-4, -4,  7, -2 * vc]], # D
     ]),

    XOR=np.array(
        [[[ 0, -1, -1,      vc],  # B
          [ 0,  1,  1,      vc],  # C
          [ 0, -1,  1,     -vc],  # E
          [ 0,  1, -1,     -vc],  # F
          [ 6,  0, -1,       0]], # D

         [[-1,  0, -1,      vc],  # B
          [ 1,  0,  1,      vc],  # C
          [-1,  0,  1,     -vc],  # E
          [ 1,  0, -1,     -vc],  # F
          [ 0,  6, -1,       0]], # D

         [[-1, -1,  0,      vc],  # B
          [ 1,  1,  0,      vc],  # C
          [-1,  1,  0,     -vc],  # E
          [ 1, -1,  0,     -vc],  # F
          [-1, -1,  7,       0]], # D
     ]),
)


class SoGate(ThreeNodeCircuitComponent):
    """
    v1--+   +--v2
        |   |
 G--DCM-+   +-DCM--G
        |   |
        R   R
        +_+_+
          |
  G--DCM--+
          |
          vo
    """
    def __init__(self, crkt, gate_params):
        self.id = get_next_id('DCM')
        super(SoGate, self).__init__(crkt)

        r_off = 2

        gnd = crkt._ground

        self.dcm1 = DynamicCorrectionModule(crkt, gate_params[0], [self.node1, self.node2, self.node3])
        self.dcm2 = DynamicCorrectionModule(crkt, gate_params[1], [self.node1, self.node2, self.node3])
        self.dcmo = DynamicCorrectionModule(crkt, gate_params[2], [self.node1, self.node2, self.node3])

        self.r1 = Resistor(crkt, r_off)
        self.r2 = Resistor(crkt, r_off)

        # hook up dcm's
        self.dcm1.positive = self.node1
        self.dcm2.positive = self.node2
        self.dcmo.positive = self.node3

        # hook up resistors
        self.r1.negative = self.node1
        self.r1.positive = self.node3
        self.r2.negative = self.node2
        self.r2.positive = self.node3

    def has_current_variable(self):
        return False

    def has_hidden_variable(self):
        return False

class SoAnd(SoGate):
    def __init__(self, crkt):
        super(SoAnd, self).__init__(crkt, GATE_PARAMS['AND'])

class SoOr(SoGate):
    def __init__(self, crkt):
        super(SoOr, self).__init__(crkt, GATE_PARAMS['OR'])

class SoXor(SoGate):
    def __init__(self, crkt):
        super(SoXor, self).__init__(crkt, GATE_PARAMS['XOR'])
