import pytest
import numpy as np

import utils as tu
import gefera as gf

coords = gf.kep.Kepler.coords
flux = gf.phot.flux
flux_ng = gf.phot.flux_ng

np.random.seed(42)

params = [
    'ab', 'tb', 'eb', 'pb', 'wb', 'ib',      
    'am', 'tm', 'em', 'pm', 'om', 'wm', 'im', 'mm'
]

def test_solver():

    n = 100
    M = np.random.rand(n) * np.pi * 2
    ecc = np.random.rand(n)
    ecc[-1] = 0.99999999
    ecc[-2] = 0.0

    tu.RPP_vs_newton(M, ecc)
    
def test_gradient():
    
    for _ in range(10):
        map(
            lambda p: tu.gradient(p), 
            params
        )
        
def test_coords():
    
    argdict = tu.random_args()
    argdict['ab'] = 0.0
    am = argdict['am']
    argdict['eb'] = 0.0
    argdict['em'] = 0.0
    
    t = np.linspace(0, 100, 10)
    xp, yp, zp, xm, ym, zm = coords(t, argdict)
    assert np.all(
        np.isclose(
            (xm - xp)**2 + (ym - yp)**2 + (zm - zp)**2, 
            am ** 2, 
            atol=1e-10
        )
    )
    
    argdict = tu.random_args()
    argdict['am'] = 0.0
    ab = argdict['ab']
    argdict['eb'] = 0.0
    argdict['em'] = 0.0
    
    t = np.linspace(0, 100, 10)
    xp, yp, zp, xm, ym, zm = coords(t, argdict)
    assert np.all(
        np.isclose(
            xm**2 + ym**2 + zm**2, 
            ab ** 2, 
            atol=1e-10
        )
    )
    assert np.all(
        np.isclose(
            xp**2 + yp**2 + zp**2, 
            ab ** 2, 
            atol=1e-10
        )
    )

def test_flux():
    
    q1, q2 = np.random.rand(2)
    u1 = 2 * np.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2 * q2)
    
    rp, rm = 0.3, 0.2
    coords = [
        [0.95, 1, 1.15, 0.75],
        [0.95, 1, 0.7, 0.75],
        [0.8, 0.8, 1.15, 0.75],
        [0.4, 0.4, 0.4, 0.4],
        [0.4, 0.4, 0.6, 0.2],
        [0.4, 0.4, 0.65, 0.65],
        [0.55, 0.55, 0.52, 0.52],
        [0.6, 0.6, 0.4, 0.4],
        [0.65, 0.65, 0.65, 0.65],
        [0.6, 0.6, 0.75, 0.75],
        [0.8, 0.8, 0.65, 0.65],
        [0.52, 0.52, 0.7, 0.7],
        [0.88, 0.88, 0.7, 0.7],
        [0.85, 0.9, 1.05, 0.5],
        [0.5, 0.55, 0.8, 0.3],
        [0.6, 0.65, 0.9, 0.4]
    ]
    for xp, yp, xm, ym in coords:
        bp, bpm, theta = tu.change_coords(xp, yp, xm, ym)    
    
        analytic = flux(
            u1, u2, rp, rm, 
            np.array([bp]), 
            np.array([bpm]), 
            np.array([np.cos(theta)]), 
            np.array([np.sin(theta)])
        )[0][0]
    
        numerical = tu.flux_numerical(
            u1, u2, rp, rm, bp, bpm, theta, 10000
        )
        assert np.isclose(numerical, analytic, atol=1e-3)