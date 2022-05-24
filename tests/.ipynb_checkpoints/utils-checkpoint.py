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
    
    a1 = uniform_draw(1, 10)
    t1 = uniform_draw(-10, 10)
    e1 = uniform_draw(0, 0.999)
    p1 = uniform_draw(1, 100)
    w1 = uniform_draw(0, np.pi*2)
    i1 = uniform_draw(0, np.pi/2)
    
    a2 = uniform_draw(1, 10)
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
    
    a1 = uniform_draw(1, 10)
    t1 = uniform_draw(-10, 10)
    e1 = uniform_draw(0, 0.999)
    p1 = uniform_draw(10, 100)
    w1 = uniform_draw(0, np.pi*2)
    i1 = uniform_draw(0, np.pi/2)
    
    a2 = uniform_draw(0.01, 1)
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