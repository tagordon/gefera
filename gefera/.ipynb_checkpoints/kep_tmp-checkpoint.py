import numpy as np
import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper_tmp.so")

__args_sat__ = [
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

__args_conf__ = [
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
        

def impacts_sat(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(3)
    )
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict, __args_sat__)
    clib.impacts_sat.restype = None
    clib.impacts_sat(
        t, 
        *args, 
        byref(ctypes.c_int(j)), bp, bpm, theta
    )
    
    return map(np.array, (bp, bpm, theta))

def impacts_conf(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(3)
    )
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict, __args_conf__)
    clib.impacts_conf.restype = None
    clib.impacts_conf(
        t, 
        *args, 
        byref(ctypes.c_int(j)), bp, bpm, theta
    )
    
    return map(np.array, (bp, bpm, theta))
    
def grad_impacts_sat(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(3)
    )
    dbp, dbpm, dtheta = tuple(
        ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14))) 
        for a in range(3)
    )
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict, __args_sat__)
    clib.grad_impacts_sat.restype = None
    clib.grad_impacts_sat(
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

def grad_impacts_conf(t, argdict):
    
    j = len(t)
    bp, bpm, theta = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(3)
    )
    dbp, dbpm, dtheta = tuple(
        ((ctypes.c_double * j) * 13).from_buffer(np.zeros((j, 13))) 
        for a in range(3)
    )
    t = byref((ctypes.c_double * j).from_buffer(t))
    
    args = validate_elements(argdict, __args_conf__)
    clib.grad_impacts_conf.restype = None
    clib.grad_impacts_conf(
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

def coords_sat(t, argdict):
    
    clib.coords_sat.restype = None
    j = len(t)
    
    xp, yp, zp, xm, ym, zm = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(6)
    )
    t = byref((ctypes.c_double * len(t)).from_buffer(t))
    
    args = validate_elements(argdict, __args_sat__)
    
    clib.coords_sat(
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

def coords_conf(t, argdict):
    
    clib.coords_conf.restype = None
    j = len(t)
    
    xp, yp, zp, xm, ym, zm = tuple(
        (ctypes.c_double * j).from_buffer(np.zeros(j)) 
        for a in range(6)
    )
    t = byref((ctypes.c_double * len(t)).from_buffer(t))
    
    args = validate_elements(argdict, __args_conf__)
    
    clib.coords_conf(
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
