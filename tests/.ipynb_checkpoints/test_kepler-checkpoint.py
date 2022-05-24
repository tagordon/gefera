import pytest
import gefera as gf

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
    
def test_impact_hrch():
    
def test_coords_conf():
    
def test_coords_hrch():
    
def test_grad_conf():
    
def test_grad_hrch():