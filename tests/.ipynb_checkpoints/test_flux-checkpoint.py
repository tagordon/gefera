import pytest
import numpy as np

import utils as tu
import gefera as gf

flux = gf.phot.flux
flux_ng = gf.phot.flux_ng

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

def test_flux():
    
    q1, q2 = np.random.rand(2)
    u1 = 2 * np.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2 * q2)
    
    rp, rm = 0.3, 0.2
    
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
            u1, u2, rp, rm, bp, bpm, theta, 500000
        )
        assert np.isclose(numerical, analytic, atol=1e-4)
        
def test_flux_ng():
    
    q1, q2 = np.random.rand(2)
    u1 = 2 * np.sqrt(q1) * q2
    u2 = np.sqrt(q1) * (1 - 2 * q2)
    
    rp, rm = 0.3, 0.2

    for xp, yp, xm, ym in coords:
        bp, bpm, theta = tu.change_coords(xp, yp, xm, ym)    
    
        analytic = flux_ng(
            u1, u2, rp, rm, 
            np.array([bp]), 
            np.array([bpm]), 
            np.array([np.cos(theta)]), 
            np.array([np.sin(theta)])
        )[0]
    
        numerical = tu.flux_numerical(
            u1, u2, rp, rm, bp, bpm, theta, 500000
        )
        assert np.isclose(numerical, analytic, atol=1e-4)
        
def test_flux_grad():
    
    u1, u2, rp, rm = 0.5, 0.3, 0.1, 0.05
    params = [u1, u2, rp, rm]
    dx = 1e-8
    
    for xp, yp, xm, ym in coords:
        bp, bpm, theta = tu.change_coords(xp, yp, xm, ym)    
    
        lc = flux(
            params[0], params[1], params[2], params[3], 
            np.array([bp]), 
            np.array([bpm]), 
            np.array([np.cos(theta)]), 
            np.array([np.sin(theta)])
        )[0]
        exact = [lc[6], lc[7], lc[1], lc[2]]
    
        fd = np.zeros_like(params)
        for i in range(len(params)):
            
            params[i] += dx
            
            flux_plus = flux_ng(
                params[0], params[1], params[2], params[3], 
                np.array([bp]), 
                np.array([bpm]), 
                np.array([np.cos(theta)]), 
                np.array([np.sin(theta)])
            )[0]
            
            params[i] -= 2 * dx
            flux_minus = flux_ng(
                params[0], params[1], params[2], params[3], 
                np.array([bp]), 
                np.array([bpm]), 
                np.array([np.cos(theta)]), 
                np.array([np.sin(theta)])
            )[0]
        
            fd[i] = (flux_plus - flux_minus) / (2 * dx)
    
        assert np.allclose(fd, exact, rtol=1e-5)