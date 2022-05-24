import numpy as np
import ctypes
from ctypes import byref

import os
import fnmatch

path, _ = os.path.split(__file__)
libname = fnmatch.filter(os.listdir(path), 'keplib*.so')[0]
clib = ctypes.CDLL(path + "/" + libname)

__args__ = [
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

def validate_elements(argdict, arglist):
    
    if set(argdict.keys()) == set(arglist):
        
        argdict = {
            k: byref(ctypes.c_double(v)) 
            for k, v in argdict.items()
        }
        idm = {v: i for i, v in enumerate(arglist)}
        args = tuple(
            dict(
                sorted(
                    argdict.items(), 
                    key=lambda p: idm[p[0]])
            ).values()
        )
        return args
        
    else:
        raise ValueError(
            "required parameters are: ", 
            arglist
        )
        
def RPP(M, ecc):
    
    j = len(M)
    M = byref((ctypes.c_double * j).from_buffer(M))
    
    for e in ecc:
        e = byref(ctypes.c_double(e))
        cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
        RPP = clib.kepler_solve_RPP(M, e, cosf, sinf, byref(ctypes.c_int(j)))
        return np.array(cosf), np.array(sinf)
    
def newton(M, ecc):
    
    j = len(M)
    M = byref((ctypes.c_double * j).from_buffer(M))
    
    for e in ecc:
        e = byref(ctypes.c_double(e))
        cosf = (ctypes.c_double * j).from_buffer(np.zeros(j))
        sinf = (ctypes.c_double * j).from_buffer(np.zeros(j))
        newton = clib.kepler_solve(M, e, cosf, sinf, byref(ctypes.c_int(j)))
        return np.array(cosf), np.array(sinf)

class Kepler:
    
    def __init__(self, args, lib):
        
        self.lib = lib
        self.args = args

    def impacts(self, t, argdict):
    
        j = len(t)
        bp, bpm, theta = tuple(
            (ctypes.c_double * j).from_buffer(np.zeros(j)) 
            for a in range(3)
        )
        t = byref((ctypes.c_double * j).from_buffer(t))
    
        args = validate_elements(argdict, self.args)
        self.lib.impacts.restype = None
        self.lib.impacts(
            t, 
            *args, 
            byref(ctypes.c_int(j)), bp, bpm, theta
        )
    
        return map(np.array, (bp, bpm, theta))
    
    def grad_impacts(self, t, argdict):
    
        j = len(t)
        bp, bpm, theta = tuple(
            (ctypes.c_double * j).from_buffer(np.zeros(j)) 
            for a in range(3)
        )
        dbp, dbpm, dtheta = tuple(
            (
                (ctypes.c_double * j) * len(self.args)
            ).from_buffer(np.zeros((j, len(self.args)))) 
            for a in range(3)
        )
        t = byref((ctypes.c_double * j).from_buffer(t))
    
        args = validate_elements(argdict, self.args)
        self.lib.grad_impacts.restype = None
        self.lib.grad_impacts(
            t, 
            *args, 
            byref(ctypes.c_int(j)), 
            bp, 
            bpm, 
            theta, 
            dbp, 
            dbpm, 
            dtheta
        )
        
        return map(
            np.array, 
            (bp, bpm, theta, dbp, dbpm, dtheta)
        )

    def coords(self, t, argdict):
    
        self.lib.coords.restype = None
        j = len(t)
    
        xp, yp, zp, xm, ym, zm = tuple(
            (ctypes.c_double * j).from_buffer(np.zeros(j)) 
            for a in range(6)
        )
        t = byref((ctypes.c_double * len(t)).from_buffer(t))
    
        args = validate_elements(argdict, self.args)
    
        self.lib.coords(
            t, 
            *args, 
            byref(ctypes.c_int(j)), 
            xp, yp, zp, 
            xm, ym, zm
        )
    
        return map(
            np.array, 
            (xp, yp, zp, xm, ym, zm)
        )
