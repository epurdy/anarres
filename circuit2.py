"""This circuit simulator handles non-linear systems with extra variables. 

This implementation is based on Najm's "Circuit Simulation", although I am
ignoring half of his advice and he should bear none of the blame for its
inadequacies.

TODO: We need some way of grouping variables or sth so that we can have numpy
handle things efficiently.
"""

from abc import ABCMeta, abstractmethod
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import scipy

import dcm
import vcdcg
from util import get_next_id


R_ON = 0.1
R_OFF = 0.2


class Var(object):
    """A variable consists of an id and a quantity. The three quantities currently
    used are 'v' for voltage (measured at nodes), 'i' for current (measured
    through some, but not all, components), and 'h' for hidden variables such
    as the memory of a memristor (measured at the component in question).
    """
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


class Stamp(object):
    """A stamp is a concept used in most circuit simulators. A stamp represents the
    contributions of a circuit element to the differential algebraic system
    governing the circuit.

    This DAS has the form

        A(x) x + B(x) x' = c(x)

    where A, B are n x n matrices and c is a n-dimensional vector.

    For some circuits containing only simple elements, A, B, and c are
    constants. In general, they depend on things like voltage, current, and
    hidden variables.
    """
    pass


class FuncStamp1(Stamp):
    """A dynamic stamp intended to populate c(x)."""
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
    """A dynamic stamp intended to populate A(x) or B(x)."""
    def __init__(self, ivar, jvar, fn, input_vars):
        self.ivar = ivar
        self.jvar = jvar
        self.fn = fn
        self.input_vars = input_vars

    def __repr__(self):
        return 'FuncStamp2({}, {}, {})'.format(self.ivar, self.jvar, self.fn.__name__)

class ConstStamp1(Stamp):
    """A constant stamp intended to populate c(x)."""
    def __init__(self, ivar, val):
        self.ivar = ivar
        self.val = val

    def __repr__(self):
        return 'ConstStamp1({}, {})'.format(self.ivar, self.val)

    def __call__(self, x, var_lut):
        return self.val

class ConstStamp2(ConstStamp1):
    """A constant stamp intended to populate A(x) or B(x)."""
    def __init__(self, ivar, jvar, val):
        self.ivar = ivar
        self.jvar = jvar
        self.val = val

    def __repr__(self):
        return 'ConstStamp2({}, {}, {})'.format(self.ivar, self.jvar, self.val)


class MnaEquation(object):
    """The differential algebraic system governing the behavior of the circuit, and
    methods for solving it.

    This DAS has the form

        A(x) x + B(x) x' = c(x)

    where A, B are n x n matrices and c is a n-dimensional vector.

    Currently, the only method for solving is the Backwards Euler method.
    """
    def __init__(self, A, B, c, rev_vdict, rev_idict, rev_hdict):
        self.Astamps = A
        self.Bstamps = B
        self.cstamps = c
        self.rev_vdict = rev_vdict
        self.rev_idict = rev_idict
        self.rev_hdict = rev_hdict
        self.nv = len(self.rev_vdict)
        self.ni = len(self.rev_idict)
        self.nh = len(self.rev_hdict)
        self.n = self.nv + self.ni + self.nh

        self.sparse = False

        # construct a lookup table for variables so that we know what
        # coordinate each corresponds to.
        self.var_lut = {}
        for k in self.rev_vdict:
            self.var_lut[Var(k, 'v')] = self.rev_vdict[k]
        for k in self.rev_idict:
            self.var_lut[Var(k, 'i')] = self.rev_idict[k]
        for k in self.rev_hdict:
            self.var_lut[Var(k, 'h')] = self.rev_hdict[k]

    def _mat(self, stamps, x):
        """Populate a matrix given a set of stamps and an x to pass in to the stamps.
        """
        if self.sparse:
            rv = scipy.sparse.dok_matrix((self.n, self.n), dtype=np.float32)
        else:
            rv = np.zeros((self.n, self.n))

        for stamp in stamps:
            i = self.var_lut[stamp.ivar]
            j = self.var_lut[stamp.jvar]

            # don't include ground
            if i is None or j is None:
                continue

            val = stamp(x, self.var_lut)
            rv[i, j] += val

        if self.sparse:
            rv = rv.tocoo()

        return rv

    def Amat(self, x):
        """Compute the matrix A"""
        return self._mat(self.Astamps, x)

    def Bmat(self, x):
        """Compute the matrix B"""
        return self._mat(self.Bstamps, x)

    def cvec(self, x):
        """Compute the vector c by populating the vector using self.cstamps."""
        rv = np.zeros(self.n)

        for stamp in self.cstamps:
            i = self.var_lut[stamp.ivar]

            # don't include ground
            if i is None:
                continue

            val = stamp(x, self.var_lut)
            rv[i] += val

        return rv

    def simulate_be(self, tmax, dt, x0=None, vars=None):
        """Simulate the equations from a given starting point using the Backwards Euler
        method.

        Args:
          tmax: Amount of time to simulate for.
          dt: The time step
          x0: Initial conditions (optional). If not provided, we will start
              from a state with zero current, very small random voltages (to
              break symmetries), and random hidden variables distributed in [0,
              1]
          vars: A list of variables whose values we want to know (optional). If
                not provided, all variables will be returned.
        """
        # set initial conditions
        if x0 is None:
            x = np.concatenate(
                [1e-6 * np.random.randn(self.nv),
                 np.zeros(self.ni),
                 np.random.random(self.nh)])
        else:
            x = np.array(x0)

        # figure out which indices of x to expose to the caller 
        if vars is None:
            indices = range(self.n)
        else:
            indices = [self.var_lut[k] for k in vars]

        # the values to return at the end
        xs = [x[indices]]

        t = 0.0
        while t < tmax:
            print t

            oldx = x
            def f(x):
                # A(x) x + B(x) x' = c(x)
                # A(x) x + B(x) * (x - oldx)/dt = c(x)
                # (A(x) + B(x)/dt) x = c(x) + B(x)/dt oldx
                A = self.Amat(x)
                B = self.Bmat(x)
                c = self.cvec(x)
                return A.dot(x) + B.dot(x - oldx) / dt - c

            # use Newton's method to find the solution to the DAS using the BE
            # discretization
            x, info, ier, mesg = scipy.optimize.fsolve(f, oldx, full_output=True)
            print np.linalg.norm(f(x))
            print mesg

            # remember the variables we wish to remember, and report their current values
            xs.append(list(x[indices]))
            print zip(vars, xs[-1])

            # step time forward
            t += dt

        xs = np.array(xs)

        return xs


class Circuit(object):
    """A circuit. This class is mainly used to build up a circuit in an
    object-oriented way. When the circuit is ready, call
    assemble_mna_equation() to get the MNA equation that will allow you to
    simulate the circuit.
    """
    def __init__(self, layout_algorithm='circular'):
        self.id = get_next_id('crkt')

        # graph where nodes are nodes and edges are two-terminal components
        self.graph = nx.MultiDiGraph()

        self.layout_algorithm = layout_algorithm

        # list of components and lookup table for components
        self.components = []
        self.component_lut = dict()

        # lookup table for nodes (necessary because nodes are constantly
        # getting collapsed with one another and stale names for nodes exist
        # outside of this class and must be handled correctly)
        self.node_lut = dict()
        
        # add ground
        self._ground = get_next_id('node')
        self.add_node(self._ground)

    def draw(self):
        # graph is significantly easier to understand without the ground node
        graph = self.graph.copy()
        graph.remove_node(self._ground)

        # select a layout algorithm
        if self.layout_algorithm == 'spring':
            pos = nx.drawing.layout.spring_layout(self.graph)
        elif self.layout_algorithm == 'circular':
            pos = nx.drawing.layout.circular_layout(self.graph)
        elif self.layout_algorithm == 'shell':
            pos = nx.drawing.layout.shell_layout(self.graph)
        elif self.layout_algorithm == 'spectral':
            pos = nx.drawing.layout.spectral_layout(graph)
        else:
            assert False, 'unrecognized layout algorithm: %s' % self.layout_algorithm

        # figure out the edge labels
        edge_labels = defaultdict(lambda: '')
        for edge in graph.edges():
            edge_data = graph.get_edge_data(*edge)
            for idx in edge_data:
                edge_labels[edge] += " " + edge_data[idx]['data']

        # draw the graph
        nx.draw_networkx(graph, pos, with_labels=True)
        nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
        
    def add_node(self, node):
        """Add a node to the circuit"""
        self.graph.add_node(node)
        self.node_lut[node] = node
        
    def add_component2(self, comp, nodea=None, nodeb=None):
        """Add a two-terminal component to the circuit"""
        if nodea is not None:
            self.add_node(nodea)

        if nodeb is not None:
            self.add_node(nodeb)

        if nodea is not None and nodeb is not None:
            self.graph.add_edge(nodea, nodeb, data=comp.id)

        self.components.append(comp)
        self.item_lut[comp.id] = comp

    def normalize(self):
        """Make sure that all node ids contained anywhere are fresh.
        """
        redo = False
        for comp in self.components:
            redo = comp.normalize() or redo

        # Sometimes we have to do multiple rounds of freshening.
        if redo:
            self.normalize()
        
    def merge(self, nodea, nodeb):
        """Merge two nodes"""

        nodea = self.node_lut[nodea]
        nodeb = self.node_lut[nodeb]

        if nodea > nodeb:
            nodea, nodeb = nodeb, nodea

        self.node_lut[nodeb] = nodea

        self.graph = nx.contracted_nodes(self.graph, nodea, nodeb, self_loops=False)

    def assemble_mna_equation(self):
        """Create the MnaEquation object needed to simulate the circuit.
        """
        self.normalize()

        # figure out the variables
        vvars = [x for x in self.graph.nodes() if x != self._ground]
        ivars = [x.id for x in self.components if x.has_current_variable()]
        hvars = [x.id for x in self.components if x.has_hidden_variable()]
        nv = len(vvars)
        ni = len(ivars)
        nh = len(hvars)
        rev_vdict = {id:num for num, id in enumerate(vvars)}
        rev_idict = {id:num + nv for num, id in enumerate(ivars)}
        rev_hdict = {id:num + nv + ni for num, id in enumerate(hvars)}

        # insert a special fake variable for ground voltage, which is always
        # defined to be zero
        rev_vdict[self._ground] = None

        # A x + B x' = c        
        Astamps = []
        Bstamps = []
        cstamps = []

        for component in self.components:
            Astamps.extend(component.Astamp())
            Bstamps.extend(component.Bstamp())
            cstamps.extend(component.cstamp())

        return MnaEquation(A=Astamps, B=Bstamps, c=cstamps, 
                           rev_vdict=rev_vdict, 
                           rev_idict=rev_idict,
                           rev_hdict=rev_hdict)

        
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
