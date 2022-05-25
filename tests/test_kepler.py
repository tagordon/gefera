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
            x2 * x2 + y2 * y2,
            atol=1e-10
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
            x2 * x2 + y2 * y2,
            atol=1e-10
        )
    )
    
def test_grad_conf():
    d = np.array(tu.random_args_conf(dictionary=False))
    t = np.linspace(0, np.max([d[3], d[9]]), 10)
    dx = 1e-8

    o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
    o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
    sys = gf.ConfocalSystem(o1, o2)

    bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
    dbp_fd = np.zeros_like(dbp)
    dbpm_fd = np.zeros_like(dbpm)
    dtheta_fd = np.zeros_like(dtheta)

    for i in range(len(d)):
        d[i] += dx
        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
        sys_plus = gf.ConfocalSystem(o1, o2)
    
        d[i] -= 2 * dx
        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
        sys_minus = gf.ConfocalSystem(o1, o2)
    
        bp_plus, bpm_plus, theta_plus = sys_plus.impacts(t)
        bp_minus, bpm_minus, theta_minus = sys_minus.impacts(t)
    
        dbp_fd[i, :] = (bp_plus - bp_minus) / (2 * dx)
        dbpm_fd[i, :] = (bpm_plus - bpm_minus) / (2 * dx)
        dtheta_fd[i, :] = (theta_plus - theta_minus) / (2 * dx)
    
        assert np.all(np.isclose(dbp[i], dbp_fd[i], atol=1e-3))
        assert np.all(np.isclose(dbpm[i], dbpm_fd[i], atol=1e-3))
        assert np.all(np.isclose(dtheta[i], dtheta_fd[i], atol=1e-3))
    
def test_grad_hrch():
    d = np.array(tu.random_args_hrch(dictionary=False))
    t = np.linspace(0, np.max([d[3], d[9]]), 10)
    dx = 1e-8

    o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
    o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
    sys = gf.HierarchicalSystem(o1, o2)

    bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
    dbp_fd = np.zeros_like(dbp)
    dbpm_fd = np.zeros_like(dbpm)
    dtheta_fd = np.zeros_like(dtheta)

    for i in range(len(d)):
        d[i] += dx
        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
        sys_plus = gf.HierarchicalSystem(o1, o2)
    
        d[i] -= 2 * dx
        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
        sys_minus = gf.HierarchicalSystem(o1, o2)
    
        bp_plus, bpm_plus, theta_plus = sys_plus.impacts(t)
        bp_minus, bpm_minus, theta_minus = sys_minus.impacts(t)
    
        dbp_fd[i, :] = (bp_plus - bp_minus) / (2 * dx)
        dbpm_fd[i, :] = (bpm_plus - bpm_minus) / (2 * dx)
        dtheta_fd[i, :] = (theta_plus - theta_minus) / (2 * dx)
    
    assert np.all(np.isclose(dbp[i], dbp_fd[i], atol=1e-3))
    assert np.all(np.isclose(dbpm[i], dbpm_fd[i], atol=1e-3))
    assert np.all(np.isclose(dtheta[i], dtheta_fd[i], atol=1e-3))