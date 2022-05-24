import numpy as np
import ctypes
from ctypes import byref
import os
import fnmatch

path, _ = os.path.split(__file__)
libname = fnmatch.filter(os.listdir(path), 'photlib*.so')[0]
clib = ctypes.CDLL(path + "/" + libname)

def flux(c1, c2, rp, rm, bp, bpm, cth, sth):
    
    j = len(bp)
    
    bp, bpm, cth, sth = map(
        lambda a: byref((ctypes.c_double * j).from_buffer(a)), 
        (bp, bpm, cth, sth, )
    )
    
    rp, rm, c1, c2 = map(
        lambda a: byref(ctypes.c_double(a)), 
        (rp, rm, c1, c2)
    )
    
    lc = ((ctypes.c_double * 8) * j).from_buffer(np.zeros((8, j)))
    j = byref(ctypes.c_int(j))
    clib.flux.restype = None
    clib.flux(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j)
    return np.array(lc)

def flux_ng(c1, c2, rp, rm, bp, bpm, cth, sth):
    
    j = len(bp)
    
    bp, bpm, cth, sth = map(
        lambda a: byref((ctypes.c_double * j).from_buffer(a)), 
        (bp, bpm, cth, sth, )
    )
    
    rp, rm, c1, c2 = map(
        lambda a: byref(ctypes.c_double(a)), 
        (rp, rm, c1, c2)
    )
    
    lc = (ctypes.c_double * j).from_buffer(np.zeros(j))
    j = byref(ctypes.c_int(j))
    clib.flux_ng.restype = None
    clib.flux_ng(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j)
    return np.array(lc)