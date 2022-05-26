import numpy as np
from copy import copy
import gefera as gf
import ctypes
from ctypes import byref

impacts = gf.kep.Kepler.impacts
grad_impacts = gf.kep.Kepler.grad_impacts
RPP = gf.kep.RPP
newton = gf.kep.newton

uniform_draw = lambda l, h: np.random.rand() * (h - l) + l

def random_args_conf(dictionary=True):
    
    keys = [
        'a1', 
        't1', 
        'e1', 
        'p1', 
        'w1', 
        'i1',      
        'a2', 
        't2', 
        'e2', 
        'p2', 
        'o2', 
        'w2', 
        'i2' 
    ]
    
    a1 = uniform_draw(1, 500)
    t1 = uniform_draw(-10, 10)
    e1 = uniform_draw(0, 0.999)
    p1 = uniform_draw(1, 100)
    w1 = uniform_draw(0, np.pi*2)
    i1 = uniform_draw(0, np.pi/2)
    
    a2 = uniform_draw(1, 500)
    t2 = uniform_draw(-10, 10)
    e2 = uniform_draw(0, 0.999)
    p2 = uniform_draw(1, 100)
    o2 = uniform_draw(0, np.pi*2)
    w2 = uniform_draw(0, np.pi*2)
    i2 = uniform_draw(0, np.pi/2)
    
    if dictionary:
        argdict = {
            k: v for k, v in zip(
                keys, 
                (a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, o2, w2, i2)
            )
        }
    
        return argdict
    
    else:
        return a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, o2, w2, i2

def random_args_hrch(dictionary=True):
    
    keys = [
        'a1', 
        't1', 
        'e1', 
        'p1', 
        'w1', 
        'i1',      
        'a2', 
        't2', 
        'e2', 
        'p2', 
        'o2', 
        'w2', 
        'i2', 
        'm2'
    ]
    
    a1 = uniform_draw(1, 500)
    t1 = uniform_draw(-10, 10)
    e1 = uniform_draw(0, 0.999)
    p1 = uniform_draw(10, 100)
    w1 = uniform_draw(0, np.pi*2)
    i1 = uniform_draw(0, np.pi/2)
    
    a2 = uniform_draw(0.1, 20)
    t2 = uniform_draw(-10, 10)
    e2 = uniform_draw(0, 0.999)
    p2 = uniform_draw(1, 10)
    o2 = uniform_draw(0, np.pi*2)
    w2 = uniform_draw(0, np.pi*2)
    i2 = uniform_draw(0, np.pi/2)
    m2 = uniform_draw(0, 1)
    
    if dictionary:
        argdict = {
            k: v for k, v in zip(
                keys, 
                (a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, o2, w2, i2, m2)
            )
        }
    
        return argdict
    
    else:
        return a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, o2, w2, i2, m2
    
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