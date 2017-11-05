"""
We use here for ground the symbol:
===
 =

The dynamic correction module circuit is as follows:

 +----VCVG_B--B--M_B--+--M_C--C--VCVG_C----+
 |                    |                    |
===                   |                   ===
 =                    |                    =
                      |
 +----VCVG_D--D---R---+--------------------------------------A------------
 |                    |               
===                   |               
 =                    |               
                      |
 +----VCVG_E--E--M_E--+--M_F--F--VCVG_F----+
 |                                         |
===                                       ===
 =                                         =

- A,B,C,D,E,F are simply nodes at which we measure the voltage
- VCVG_X is a voltage-controlled voltage generator
- M_X is a memristor
- R is a standard resistor
"""

import numpy as np
from vcdcg import smooth_step_function_np

def memristor_conductance(x, r_on=1.0, r_off=2.0):
    resistance = r_on + (r_off - r_on) * x
    assert (resistance > 0).all()
    return 1.0 / resistance


def derivative_of_memristor_voltage(i_c, C=0.001):
    assert C > 0
    return i_c / C
    

def memristor_coupling(x, v_m, k=1000.0):
    """h(x, v_m)"""
    x = np.clip(x, 0.0, 1.0)
    term1 = (1 - np.exp(-k * x)) * smooth_step_function_np(v_m)
    term2 = (1 - np.exp(-k * (1 - x))) * smooth_step_function_np(-v_m)
    return term1 + term2


def derivative_of_memristor_variables(v_m, x, alpha=1.0):
    return -alpha * memristor_coupling(x, v_m) * memristor_conductance(x) * v_m


def test_single_dcm():
    """Test a partial AND gate where we provide the input voltage v2 and the output
    voltage vo, and v1 has both a vcdcg and a dcm hooked up to it."""

    # run nthings of these AND gates in parallel
    nthings = 10

    # the real variables
    v_m = 20 * np.random.random(nthings, 5) - 10
    x = np.random.random(nthings, 5)
    i_dcg = 40 * np.random.random(nthings) - 20
    s = np.random.random(nthings)

    v_datapoints = []

    R = 1.0
    h = 0.001
    for i in xrange(10 * 1000):
        print i

        # compute the aux variables
        i_c = i_tot - v_m / memristor_conductance(x)
        v_dcg = v_dcg_given_v_m(v_m)
        v_datapoints.append(v_dcg)

        # compute derivatives
        d_v_m = derivative_of_memristor_voltage(i_c)
        d_x = derivative_of_memristor_variables(v_m, x)
        d_i_dcg = derivative_of_dcg_currents(v_dcg, i_dcg, s, delta_s=1e-3, gamma=1.0, v_c=10.0)
        d_s = derivative_of_dcg_variables(i_dcg, s, i_min=5.0, i_max=15.0, 
                                          k_s=1.0, k_i=10.0, delta_i=0.001)

        # integrate forward using the dumbest possible method
        v_m += h * d_v_m
        x += h * dx
        i_dcg += h * d_i_dcg
        s += h * d_s
        s_datapoints.append(s)

        print np.any(np.abs(i_dcg) > 0.1), np.all(np.abs(i_dcg) < 19.9)
        #print s

    plt.plot(v_datapoints)
    plt.savefig('dcm.png')


test_single_dcm()
