import numpy as np
from timeit import default_timer as timer
import ctypes
import os
import fnmatch

__all__ = ['PrimaryOrbit', 'SatelliteOrbit', 'ConfocalOrbit']

class Orbit:
    
    """Parent class of all orbits. This class does not contain
    all the necessary attributes for all orbits. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    
    def __init__(self, a, t, e, p, w, i):
        
        self.a = a
        self.t = t
        self.e = e
        self.p = p
        self.w = w
        self.i = i
        
    def pdict(self):
        
        return vars(self)
    
class PrimaryOrbit(Orbit):
    
    """
    A heliocentric orbit for the primary body in the system. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    
    def __init__(self, a, t, e, p, w, i):
        super().__init__(a, t, e, p, w, i)
        
    def pdict(self):
        
        return {k + '1': v for k, v in vars(self).items()}
    
class SatelliteOrbit(Orbit):
    
    """
    The orbit of the moon around the planet.
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        o: Longitude of ascending node (in radians)
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
        m: Moon/planet mass ratio
    """
    
    def __init__(self, a, t, e, p, o, w, i, m):
        super().__init__(a, t, e, p, w, i)
        self.o = o
        self.m = m
        
    def pdict(self):
        
        return {k + '2': v for k, v in vars(self).items()}
        
class ConfocalOrbit(Orbit):
    
    """
    A second heliocentric orbit. 
    
    Args:
        a: Semimajor axis
        t: Time of periastron passage
        e: Eccentricity
        p: Period
        o: Longitude of ascending node (in radians)
        w: Argument of periastron (in radians)
        i: Inclination (in radians)
    """
    
    def __init__(self, a, t, e, p, o, w, i):
        super().__init__(a, t, e, p, w, i)
        self.o = o
        
    def pdict(self):
        
        return {k + '2': v for k, v in vars(self).items()}