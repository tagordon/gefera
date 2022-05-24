import pytest
import numpy as np
import gefera as gf

import utils as tu

RPP = gf.kep.RPP
newton = gf.kep.newton

def test_solver():

    n = 100
    M = np.random.rand(n) * np.pi * 2
    ecc = np.random.rand(n)
    ecc[-1] = 0.99999999
    ecc[-2] = 0.0
    
    cosf_RPP, sinf_RPP = RPP(M, ecc)
    cosf_newton, sinf_newton = newton(M, ecc)
    assert np.all(np.isclose(np.array(sinf_RPP), np.array(sinf_newton), atol=1e-12))
    assert np.all(np.isclose(np.array(cosf_RPP), np.array(cosf_newton), atol=1e-12))

def test_impact_conf():
    
    d = tu.random_args_conf()
    o1 = gf.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.ConfocalOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'])
    sys = gf.ConfocalSystem(o1, o2)
    
    t = np.linspace(0, np.max([d['p1'], d['p2']]), 1000)
    (x1, y1, z1), (x2, y2, z2) = sys.coords(t)
    bp, bpm, theta = sys.impacts(t)
    
    # check consistency with law of cosines
    assert np.all(
        np.isclose(
            bpm * bpm + bp * bp - 2 * bp * bpm * np.cos(theta),
            x2 * x2 + y2 * y2
        )
    )
    
def test_impact_hrch():
    
    d = tu.random_args_hrch()
    o1 = gf.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.SatelliteOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'], d['m2'])
    sys = gf.HierarchicalSystem(o1, o2)
    
    t = np.linspace(0, np.max([d['p1'], d['p2']]), 1000)
    (x1, y1, z1), (x2, y2, z2) = sys.coords(t)
    bp, bpm, theta = sys.impacts(t)
    
    # check consistency with law of cosines
    assert np.all(
        np.isclose(
            bpm * bpm + bp * bp - 2 * bp * bpm * np.cos(theta),
            x2 * x2 + y2 * y2
        )
    )
    
def test_coords_conf():
    pass
    
def test_coords_hrch():
    pass
    
def test_grad_conf():
    pass
    
def test_grad_hrch():
    pass