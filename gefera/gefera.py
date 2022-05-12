import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from timeit import default_timer as timer
import ctypes

from .kep import Kepler
from .phot import flux, flux_ng

__all__ = ['PrimaryOrbit', 'SatelliteOrbit', 'ConfocalOrbit']

au_r = 215.03215567054764

conflib = ctypes.CDLL("../builddir/conf.so")
hrchlib = ctypes.CDLL("../builddir/hrch.so")

confargs = [
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

hrchargs = [
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
        
class System:
    
    """Base class of all systems. This class should not 
    be instantiated directly. Use ConfocalSystem or 
    HierarchicalSystem instead.
    
    Args:
        o1 (Orbit): The orbit of the first body 
        o2 (Orbit): The orbit of the second body 
    """
    
    def __init__(self, o1, o2):
        
        self.o1 = o1
        self.o2 = o2
        self.pdict = {**o1.pdict(), **o2.pdict()}
            
    def coords(self, t):
        
        """
        Get the coordinates of the planet and moon.
        
        Args:
            t: Times at which the coordinates should be computed.
            
        Returns:
            pc: Coorinates of the planet as a tuple of arrays (x, y, z)
            mc: Coordinates of the moon as a tuple of arrays (x, y, z)
        """
        
        x1, y1, z1, x2, y2, z2 = self.kep.coords(t, self.pdict)
        return ((x1 * au_r, y1 * au_r, z1 * au_r), 
                (x2 * au_r, y2 * au_r, z2 * au_r))
    
    def impacts(self, t, grad=False):
        
        if grad:
            bp, bpm, theta, dbp, dbpm, dtheta = self.kep.grad_impacts(
                t, 
                self.pdict
            )
            return bp, bpm, theta, dbp, dbpm, dtheta
        else:
            bp, bpm, theta = self.kep.impacts(t, self.pdict)
            return bp, bpm, theta
        
    def phot(self, t, u1, u2, r1, r2, bp, bpm, theta, grad=False):
        
        if grad:
            lc = flux(
                u1, 
                u2, 
                r1, 
                r2, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            ).T
        else:
            lc = flux_ng(
                u1, 
                u2, 
                r1, 
                r2, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            )
        
        return lc
    
    def lightcurve(self, t, u1, u2, r1, r2, grad=False, integrate=None, dt=None):
        
        """
        Get the lightcurve resulting from a transit of the moon/planet system.
        
        Args: 
            t: Times at which the flux should be computed
            r1: Radius of the body in the PrimaryOrbit
            r2: Radius of the body in the SecondaryOrbit or ConfocalOrbit
            u1: The first limb-darkening parameter
            u2: The second limb-darkening parameter
            grad (bool): If True, compute the gradient of the lightcurve.
                Default is False. 
                
        Returns:
            lc: The lightcurve
            grad: A dict containing the derivatives of the 
                lightcurve with respect to each of the input parameters.
            
        """
        if integrate == 'simpson':
            ta = t - dt / 2
            tb = t + dt / 2
            
            if grad:
                f, g = self.lightcurve(t, u1, u2, r1, r2, grad=True)
                fa, ga = self.lightcurve(ta, u1, u2, r1, r2, grad=True)
                fb, gb = self.lightcurve(tb, u1, u2, r1, r2, grad=True)
                lc = (fa + 4 * f + fb) / 6
                grad = {
                        ka: (va + 4 * v + vb) / 6
                        for (k, v), (ka, va), (kb, vb) in zip(g.items(), ga.items(), gb.items())
                }
                return lc, grad
            else:
                f = self.lightcurve(t, u1, u2, r1, r2, grad=False)
                fa = self.lightcurve(ta, u1, u2, r1, r2, grad=False)
                fb = self.lightcurve(tb, u1, u2, r1, r2, grad=False)
                return (fa + 4 * f + fb) / 6
        
        if integrate == 'trapezoid':
            
            if dt is None:
                dt = np.diff(t)
                tt = np.zeros(len(t) + 1)
                tt[1:-1] = t[1:] - dt / 2
                tt[0] = t[0] - dt[0] / 2
                tt[-1] = t[-1] + dt[-1] / 2
                
                if grad:
                    f, g = self.lightcurve(tt, u1, u2, r1, r2, grad=True)
                    lc = (f[1:] + f[:-1]) / 2
                    grad = {
                        k: (v[1:] + v[:-1]) / 2
                        for k, v in g.items()
                    }
                    return lc, grad
                else:
                    f = self.lightcurve(tt, u1, u2, r1, r2, grad=False)
                    return (f[1:] + f[:-1]) / 2
            else:
            
                ta = t - dt / 2
                tb = t + dt / 2
                if grad:
                    fa, ga = self.lightcurve(ta, u1, u2, r1, r2, grad=True)
                    fb, gb = self.lightcurve(tb, u1, u2, r1, r2, grad=True)
                    lc = (fa + fb) / 2
                    grad = {
                        ka: (va + vb) / 2
                        for (ka, va), (kb, vb) in zip(ga.items(), gb.items())
                    }
                    return lc, grad
                else:
                    fa = self.lightcurve(ta, u1, u2, r1, r2)
                    fb = self.lightcurve(tb, u1, u2, r1, r2)
                    return (fa + fb) / 2
        
        if grad:
            bp, bpm, theta, dbp, dbpm, dtheta = self.kep.grad_impacts(
                t, 
                self.pdict
            )
            f = flux(
                u1, 
                u2, 
                r1, 
                r2, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            ).T
                
            lc = f[0]
            f_bp = f[3]
            f_bpm = f[4]
            f_theta = f[5]
                
            df = (
                f_bp * dbp 
                + f_bpm * dbpm 
                + f_theta * dtheta
            )
                
            grad = {
                self.argnames[i]: df[i] 
                for i in range(np.shape(df)[0])
            }
            
            # order?
            grad['r1'] = f[1]
            grad['r2'] = f[2]
            grad['u1'] = f[6]
            grad['u2'] = f[7]
            
            return lc, grad
        
        else:
            bp, bpm, theta = self.kep.impacts(
                t, 
                self.pdict
            )
            return flux_ng(
                u1, 
                u2, 
                r1, 
                r2, 
                bp, 
                bpm, 
                np.cos(theta), 
                np.sin(theta)
            )
                    
    def loglike(self, y, t, u1, u2, r1, r2, sigma, integrate=None, dt=None, grad=False, sign=1):
        
        """
        Get the log-likelihood of the lightcurve.
        
        Args:
            y: A vector of observations to compute the likelihood with 
                respect to. 
            t: Times at which the flux should be computed
            r1: Radius of the body in the PrimaryOrbit
            r2: Radius of the body in the SecondaryOrbit or ConfocalOrbit
            u1: The first limb-darkening parameter
            u2: The second limb-darkening parameter
            sigma: The standard deviation of the model
        """
        
        if grad:
            mu, jac = self.lightcurve(t, u1, u2, r1, r2, integrate=integrate, dt=dt, grad=True)
            s2 = sigma * sigma
            jac = np.array(list(jac.values()))
            ll = -0.5 * np.sum((y - mu) ** 2 / s2 + np.log(s2))
            dldsig = np.sum((y - mu) ** 2 / (s2 * sigma) - 1 / sigma)
            dldx = np.sum((y - mu) * jac / s2, axis=1)
            return sign * ll, sign * np.hstack([dldsig, dldx])
        else:
            mu = self.lightcurve(t, u1, u2, r1, r2, integrate=None, dt=None)
            s2 = sigma * sigma
            return sign * -0.5 * np.sum((y - mu) ** 2 / s2 + np.log(s2))

    def time(self, t, u1, u2, r1, r2, phot_only=False, grad=False, ntimes=1):
        
        if phot_only:
            bp, bpm, theta = self.kep.impacts(t)
            if grad:
                start = timer()
                for _ in range(ntimes):
                    lc = self.phot(t, u1, u2, r1, r2, bp, bpm, theta, grad=True)
                end = timer()
            else:
                start = timer()
                for _ in range(ntimes):
                    lc = self.phot(t, u1, u2, r1, r2, bp, bpm, theta)
                end = timer()
        else:
            if grad:
                start = timer()
                for _ in range(ntimes):
                    bp, bpm, theta, dbp, dbpm, dtheta = self.kep.impacts(t, grad=True)
                    lc = self.phot(t, u1, u2, r1, r2, bp, bpm, theta, grad=True)
                    df = (
                        lc[3] * dbp 
                        + lc[4] * dbpm 
                        + lc[5] * dtheta
                    )
                end = timer()
            else:
                start = timer()
                for _ in range(ntimes):
                    bp, bpm, theta = self.kep.impacts(t)
                    lc = self.phot(t, u1, u2, r1, r2, bp, bpm, theta)
                end = timer()
        return (end - start) / ntimes
    
class ConfocalSystem(System):
    
    """
    Represents a system with two bodies 
    orbiting a central star. 
    
    Args:
        o1 (PrimaryOrbit): The orbit of the first body 
        o2 (ConfocalOrbit): The orbit of the second body
    """
    
    def __init__(self, o1, o2):
        super().__init__(o1, o2)
        
        if not isinstance(o1, PrimaryOrbit):
            msg = ("o1 should be a PrimaryOrbit")
            raise AttributeError(msg)
        if not isinstance(o2, ConfocalOrbit):
            msg = ("o2 should be a ConfocalOrbit")
            raise AttributeError(msg)
        
        self.argnames = confargs
        self.kep = Kepler(self.argnames, conflib)
        
class HierarchicalSystem(System):
    
    """
    Represents a system with a primary body orbiting 
    a central star and a satellite orbiting the 
    primary body.
    
    Args:
        o1 (PrimaryOrbit): The orbit of the primary 
                           body around the central star
        o2 (SatelliteOrbit): The orbit of the satellite around 
                             the primary body.
    """
    
    def __init__(self, o1, o2):
        super().__init__(o1, o2)
        
        if not isinstance(o1, PrimaryOrbit):
            msg = ("o1 should be a PrimaryOrbit")
            raise AttributeError(msg)
        if not isinstance(o2, SatelliteOrbit):
            msg = ("o2 should be a SatelliteOrbit")
            raise AttributeError(msg)
        
        self.argnames = hrchargs
        self.kep = Kepler(self.argnames, hrchlib) 
        
