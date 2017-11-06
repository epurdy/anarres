#from typing import Any
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def smooth_step_function_np(x):
    x = np.minimum(x, 1.0)
    x = np.maximum(x, 0.0)
    return -20.0 * x ** 7 + 70 * x ** 6 - 84 * x ** 5 + 35 * x ** 4


def test_smooth_step_function():
    x = np.linspace(-1, 2, 301)
    plt.plot(x, smooth_step_function_np(x))
    plt.savefig('step.png')


def derivative_of_dcg_variables(i_dcg, s, 
                                i_min, i_max, 
                                k_s, k_i, 
                                delta_i):
    """f_s(i_dcg, s)"""

    term1 = -k_s * s * (s - 1) * (2 * s - 1)
    
    thing1 = (i_min * i_min - i_dcg * i_dcg) / delta_i
    thing2 = (i_max * i_max - i_dcg * i_dcg) / delta_i
    
    term2 = k_i * (1 - smooth_step_function_np(thing1).prod() 
                   - smooth_step_function_np(thing2).prod())

    # experiment
    term2 = -term2

    # print term1, term2

    return term1 + term2


def test_derivative_of_dcg_variables():
    i_dcg = np.ones(301)
    s = np.linspace(-0.2, 1.2, 301)
    deriv = derivative_of_dcg_variables(i_dcg, s, i_min=0.1, i_max=19.9, k_s=1.0, k_i=10.0, delta_i=0.001)
    plt.plot(s, deriv)
    plt.savefig('dcg_vars.png')


def bistable_fn(v, v_c):
    factor = 1.0 / v_c
    right_hump = -1 + 2 * smooth_step_function_np(factor * (v - v_c) + 0.5)
    middle_hump = 1 - 2 * smooth_step_function_np(factor * v + 0.5)
    left_hump = -1 + 2 * smooth_step_function_np(factor * (v + v_c) + 0.5)
    return np.minimum(left_hump, np.maximum(middle_hump, right_hump))

def bistable_fn2(v, v_c):
    factor = 4.0
    right_hump = np.minimum(np.maximum(factor * (v - v_c), -v_c), v_c)
    middle_hump = np.minimum(np.maximum(-factor * v, -v_c), v_c)
    left_hump = np.minimum(np.maximum(factor * (v + v_c), -v_c), v_c)
    return np.minimum(left_hump, np.maximum(middle_hump, right_hump))


def test_bistable_fn():
    vs = np.linspace(-20, 20, 10001)
    f_dcg = bistable_fn2(vs, 10)
    plt.plot(vs, f_dcg)
    plt.savefig('bistable.png')


def derivative_of_dcg_currents(v_dcg, i_dcg, 
                               s, delta_s,
                               gamma, v_c):
    rho_of_s = smooth_step_function_np((s - 0.5) / delta_s)
    rho_of_oms = smooth_step_function_np(((1 - s) - 0.5) / delta_s)

    f_dcg = -bistable_fn2(v_dcg, v_c)

    # rv = rho_of_s * f_dcg - gamma * rho_of_oms * i_dcg
    rv = rho_of_oms * f_dcg - gamma * rho_of_s * i_dcg

    # print 's', s
    # print 'v_dcg', v_dcg
    # print 'i_dcg', i_dcg
    # print 'f_dcg', f_dcg
    # print 'rv', rv

    return rv


def test_simple_system():
    """See what happens when we connect a bunch of VCDCG's in parallel. Each VCDCG
    has a 1 Ohm resistor in series with it."""
    nthings = 10

    i_dcg = 40 * np.random.random(nthings) - 20
    #i_dcg = 0.1 * (np.random.random(nthings) - 0.5)
    s = np.random.random(nthings)

    v_datapoints = []
    s_datapoints = []

    R = 1.0
    h = 0.001
    for i in xrange(10 * 1000):
        print i
        v_dcg = i_dcg * R
        v_datapoints.append(v_dcg)
        d_i_dcg = derivative_of_dcg_currents(v_dcg, i_dcg, s, delta_s=1e-3, gamma=1.0, v_c=10.0)
        d_s = derivative_of_dcg_variables(i_dcg, s, i_min=5.0, i_max=15.0, 
                                          k_s=1.0, k_i=10.0, delta_i=0.001)
        i_dcg += h * d_i_dcg
        s += h * d_s
        s_datapoints.append(s)

        print np.any(np.abs(i_dcg) > 0.1), np.all(np.abs(i_dcg) < 19.9)
        #print s

    plt.plot(v_datapoints)
    plt.savefig('dynamics.png')

    plt.figure()
    plt.plot(s_datapoints)
    plt.savefig('dynamics2.png')


if __name__ == '__main__':
    test_simple_system()
    #test_derivative_of_dcg_variables()
