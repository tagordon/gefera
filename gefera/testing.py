import numpy as np
from copy import copy
import sys
sys.path.append('../gefera')
from kep import impacts, grad_impacts
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

clib.kepler_solve.restype = None
clib.kepler_solve_RPP.restype = None

uniform_draw = lambda l, h: np.random.rand() * (h - l) + l

def random_args():
    
    keys = [
        'ab', 
        'tb', 
        'eb', 
        'pb', 
        'wb', 
        'ib',      
        'am', 
        'tm', 
        'em', 
        'pm', 
        'om', 
        'wm', 
        'im', 
        'mm'
    ]
    
    ab = uniform_draw(1, 10)
    tb = uniform_draw(-10, 10)
    eb = uniform_draw(0, 0.999)
    pb = uniform_draw(10, 100)
    wb = uniform_draw(0, np.pi*2)
    ib = uniform_draw(0, np.pi/2)
    
    am = uniform_draw(0.01, 1)
    tm = uniform_draw(-10, 10)
    em = uniform_draw(0, 0.999)
    pm = uniform_draw(1, 10)
    om = uniform_draw(0, np.pi*2)
    wm = uniform_draw(0, np.pi*2)
    im = uniform_draw(0, np.pi/2)
    mm = uniform_draw(0, 1)
    
    argdict = {
        k: v for k, v in zip(
            keys, 
            (ab, tb, eb, pb, wb, ib, am, tm, em, pm, om, wm, im, mm)
        )
    }
    
    return argdict

def gradient(p):
    d = 0.000001
    t = np.linspace(0, 100, 10)
    argdict = random_args()
    ind = np.where(np.array(list(argdict.keys())) == p)[0][0]
    pargdict = copy(argdict)
    margdict = copy(argdict)
    pargdict[p] = argdict[p] + d
    margdict[p] = argdict[p] - d
    pbp, pbpm, ptheta = impacts(t, pargdict)
    mbp, mbpm, mtheta = impacts(t, margdict)
    _, _, _, dbp, dbpm, dtheta = grad_impacts(t, argdict)
    fdiff_bp = (pbp - mbp) / (2 * d)
    fdiff_bpm = (pbpm - mbpm) / (2 * d)
    fdiff_theta = (ptheta - mtheta) / (2 * d)
    
    tol = 1e-3
    
    if p == 'pb':
        tol = 0.1
    if p == 'pm':
        tol = 0.5
    
    assert np.all(np.isclose(fdiff_bp, dbp[ind], atol=tol))
    assert np.all(np.isclose(fdiff_bpm, dbpm[ind], atol=tol))
    assert np.all(np.isclose(fdiff_theta, dtheta[ind], atol=tol))

def RPP_vs_newton(M, ecc):
    
    j = len(M)
    M = byref((ctypes.c_double * j).from_buffer(M))
    
    for e in ecc:
        e = byref(ctypes.c_double(e))
        cosf_newton = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf_newton = (ctypes.c_double * j).from_buffer(np.zeros(j))
        cosf_RPP = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf_RPP = (ctypes.c_double * j).from_buffer(np.zeros(j))
        newton = clib.kepler_solve(M, e, cosf_newton, sinf_newton, byref(ctypes.c_int(j)))
        RPP = clib.kepler_solve_RPP(M, e, cosf_RPP, sinf_RPP, byref(ctypes.c_int(j)))
        assert np.all(np.isclose(np.array(sinf_RPP), np.array(sinf_newton), atol=1e-12))
        assert np.all(np.isclose(np.array(cosf_RPP), np.array(cosf_newton), atol=1e-12))

def flux_numerical(u1, u2, rp, rm, bp, bpm, theta, nn):
    
    xp, yp = -bp, 0.0
    xm, ym = -bp + bpm * np.cos(theta), bpm * np.sin(theta)
    
    ld = lambda p: [
        1 - u1 * (1 - np.sqrt(1 - p[0]**2 - p[1]**2)) 
        - u2 * (1 - np.sqrt(1 - p[0]**2 - p[1]**2))**2 
        for p in p
    ]
    ld1 = 1 - u1 - 2 * u2
    ld2 = u1 + 2 * u2
    ld3 = u2
    
    norm = ld1 * np.pi + ld2 * 2 * np.pi / 3 + ld3 * np.pi / 2
    rp2 = rp * rp
    rm2 = rm * rm
    bpm = np.sqrt((xp - xm)**2 + (yp - ym)**2)
    
    res = nn ** 2

    if bpm < (rp + rm):
        
        if bpm + rm > rp:
            scale = 2 * (bpm + rm) + rm / 10
        else: 
            scale = 2 * rp + rm / 10
        
        n = int(np.sqrt(res * (4 * rp) ** 2))
        rat = np.sqrt(3) / 2
        nx = int(np.sqrt(n / rat))
        ny = n // nx
        xv, yv = np.meshgrid(np.arange(nx), np.arange(ny))
        xv = xv * rat
        xv[::2, :] += rat / 2
        xv = xv / len(xv) * scale + xp - scale / 2
        yv = yv / len(yv) * scale + yp - scale / 2
        points = np.concatenate((xv.flatten()[:, None], yv.flatten()[:, None]), axis=1)
    else:
        lxp = xp - rp
        uxp = xp + rp
        lxm = xm - rm
        uxm = xm + rm
    
        lyp = yp - rp
        uyp = yp + rp
        lym = ym - rm
        uym = ym + rm
    
        nump = int(np.sqrt(res * (2 * rp) ** 2))
        rat = np.sqrt(3) / 2
        nxp = int(np.sqrt(nump / rat))
        nyp = nump // nxp
        xvp, yvp = np.meshgrid(np.arange(nxp), np.arange(nyp))
        xvp = xvp * rat
        xvp[::2, :] += rat / 2
        xvp = xvp / len(xvp) * 2 * rp + xp - rp
        yvp = yvp / len(yvp) * 2 * rp + yp - rp
        pointsp = np.concatenate((xvp.flatten()[:, None], yvp.flatten()[:, None]), axis=1)
        
        numm = int(np.sqrt(res * (2 * rm) ** 2))
        rat = np.sqrt(3) / 2
        nxm = int(np.sqrt(numm / rat))
        nym = numm // nxm
        xvm, yvm = np.meshgrid(np.arange(nxm), np.arange(nym))
        xvm = xvm * rat
        xvm[::2, :] += rat / 2
        xvm = xvm / len(xvm) * 2 * rm + xm - rm
        yvm = yvm / len(yvm) * 2 * rm + ym - rm
        pointsm = np.concatenate((xvm.flatten()[:, None], yvm.flatten()[:, None]), axis=1)
        points = np.concatenate((pointsp, pointsm))
    
    dx = points[1, 0] - points[0, 0]
    dd =  dx * dx * 2 / np.sqrt(3)
    dp2 = (points[:, 0] - xp)**2 + (points[:, 1] - yp)**2
    dm2 = (points[:, 0] - xm)**2 + (points[:, 1] - ym)**2
    ds2 = points[:, 0]**2 + points[:, 1]**2
    
    occulted = ((dp2 < rp2) | (dm2 < rm2)) & (ds2 < 1)
    return - np.sum(ld(points[occulted])) * dd / norm

def change_coords(xp, yp, xm, ym):
    
    bm2 = xm**2 + ym**2
    bm = np.sqrt(bm2)
    
    bp2 = xp**2 + yp**2
    bp = np.sqrt(bp2)
    
    bpm2 = (xm - xp)**2 + (ym - yp)**2
    bpm = np.sqrt(bpm2)

    a = bp 
    b = bpm
    c = bm
        
    if a > b:
        tmp = b
        b = a
        a = tmp
    if b > c:
        mu = c - (a - b)
    else:
        mu = b - (a - c)
    
    if bpm == 0 or (a - b + c) < 1e-15:
        theta = 0.0
    elif (a - c + b) < 1e-15:
        theta = np.pi
    else:
        theta = 2 * np.arctan(np.sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
    return bp, bpm, theta
