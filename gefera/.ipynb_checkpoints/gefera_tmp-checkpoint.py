import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from timeit import default_timer as timer

#from kep_tmp import impacts, coords, grad_impacts
import kep_tmp as kep
from phot import flux, flux_ng

au_r = 215.03215567054764

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
        
        x1, y1, z1, x2, y2, z2 = self.coordfunc(t, self.pdict)
        return ((x1 * au_r, y1 * au_r, z1 * au_r), 
                (x2 * au_r, y2 * au_r, z2 * au_r))
    
    def impacts(self, t, grad=False):
        
        if grad:
            bp, bpm, theta, dbp, dbpm, dtheta = self.gimpfunc(
                t, 
                self.pdict
            )
            return bp, bpm, theta, dbp, dbpm, dtheta
        else:
            bp, bpm, theta = self.impfunc(t, self.pdict)
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
    
    def lightcurve(self, t, u1, u2, r1, r2, grad=False):
        
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
        
        if grad:
            bp, bpm, theta, dbp, dbpm, dtheta = self.gimpfunc(
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
            grad['r2'] = f[1]
            grad['r1'] = f[2]
            grad['u1'] = f[6]
            grad['u2'] = f[7]
            
            return lc, grad
        
        else:
            bp, bpm, theta = self.impfunc(
                t, 
                self.pdict
            )
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
        
    def loglike(self, y, t, u1, u2, r1, r2, sigma):
        
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
        
        mu = self.lightcurve(t, u1, u2, r1, r2)
        s2 = sigma * sigma
        return -0.5 * np.sum((y - mu) ** 2 / s2 + np.log(s2))

    def time(self, t, u1, u2, r1, r2, phot_only=False, grad=False, ntimes=1):
        
        if phot_only:
            bp, bpm, theta = self.impacts(t)
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
                    bp, bpm, theta, dbp, dbpm, dtheta = self.impacts(t, grad=True)
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
                    bp, bpm, theta = self.impacts(t)
                    lc = self.phot(t, u1, u2, r1, r2, bp, bpm, theta)
                end = timer()
        return (end - start) / ntimes
    
    def draw(self, ax, t, r1, r2, ld_params=None, cmap=plt.cm.copper, fill=True):
        
        if isinstance(t, np.ndarray):
            raise Exception("Argument t should be a scalar, not an array.")
            
        x1, y1, z1, x2, y2, z2 = self.coordfunc(np.array([t]), self.pdict)            
        
        b1 = plt.Circle(
            (x1, y1),
            radius=r1,
            color='k',
            fill=fill
        )
        b2 = plt.Circle(
            (x2, y2),
            radius=r2,
            color='k',
            fill=fill
        )
        if ld_params is None:
            star = plt.Circle((0, 0), radius=1, color=cmap(1.0), fill=True)
            ax.add_patch(star)
        
        else:
            u1, u2 = ld_params
            
            d = 0.01
            x = np.arange(-1.0, 1.0, d)
            y = np.arange(-1.0, 1.0, d)
            x, y = np.meshgrid(x, y)
            shape = np.shape(x)
            x = x.flatten()
            y = y.flatten()
            z = np.zeros(len(x))
            with np.errstate(invalid='ignore'):
                for i, (x, y) in enumerate(zip(x, y)):
                    z[i] = 1 - u1 * (1 - np.sqrt(1 - x**2 - y**2)) - u2 * (1 - np.sqrt(1 - x**2 - y**2))**2
            z = z.reshape(shape)
            
            cmap = cmap
            cmap.set_bad('white')

            edge = plt.Circle(
                (0.005, -0.0), 
                radius=1.0, 
                color='w',
                fill=False, 
                linewidth=4
            )
            
            ax.imshow(z, interpolation='bilinear', extent=(-1, 1, -1, 1), cmap=cmap)
            ax.add_patch(edge)
            
        ax.add_patch(b1)
        ax.add_patch(b2)
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_axis_off()
        
    def animate(self, fig, t, r1, r2, duration=5, ld_params=None, cmap=plt.cm.copper):
        
        t = np.ascontiguousarray(t)
        ax = fig.gca()
            
        x1, y1, z1, x2, y2, z2 = self.coordfunc(t, self.pdict)
        
        if ld_params is None:
            star = plt.Circle((0, 0), radius=1, color=cmap(1.0), fill=True)
            ax.add_patch(star)
        
        else:
            u1, u2 = ld_params
            
            d = 0.01
            x = np.arange(-1.0, 1.0, d)
            y = np.arange(-1.0, 1.0, d)
            x, y = np.meshgrid(x, y)
            shape = np.shape(x)
            x = x.flatten()
            y = y.flatten()
            z = np.zeros(len(x))
            with np.errstate(invalid='ignore'):
                for i, (x, y) in enumerate(zip(x, y)):
                    z[i] = 1 - u1 * (1 - np.sqrt(1 - x**2 - y**2)) - u2 * (1 - np.sqrt(1 - x**2 - y**2))**2
            z = z.reshape(shape)
            
            cmap = cmap
            cmap.set_bad('white')

            edge = plt.Circle(
                (0.005, -0.0), 
                radius=1.0, 
                color='w',
                fill=False, 
                linewidth=4
            )
        
            ax.imshow(z, interpolation='bilinear', extent=(-1, 1, -1, 1), cmap=cmap)
            ax.add_patch(edge)
        
        interval = int(duration * 1000 / len(t))
        
        p_patch = ax.add_patch(plt.Circle((0, 0), 0, color='k'))
        m_patch = ax.add_patch(plt.Circle((0, 0), 0, color='k'))
        
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        
        ax.set_axis_off()
        
        def init():

            p_patch.set_center((0, 0))
            p_patch.set_radius(r1)
            m_patch.set_center((0, 0))
            m_patch.set_radius(r2)
            return p_patch, m_patch
        
        def update(i):

            p_patch.set_center((x1[i], y1[i]))
            m_patch.set_center((x2[i], y2[i]))
            return p_patch, m_patch
        
        return animation.FuncAnimation(
            fig, 
            update, 
            frames=np.arange(len(t)), 
            init_func=init, 
            blit=True, 
            interval=interval
        )
    
class ConfocalSystem(System):
    
    def __init__(self, o1, o2):
        super().__init__(o1, o2)
        
        if not isinstance(o1, PrimaryOrbit):
            msg = ("o1 should be a PrimaryOrbit")
            raise AttributeError(msg)
        if not isinstance(o2, ConfocalOrbit):
            msg = ("o2 should be a ConfocalOrbit")
            raise AttributeError(msg)
            
        self.coordfunc = kep.coords_conf
        self.impfunc = kep.impacts_conf
        self.gimpfunc = kep.grad_impacts_conf
        
        self.argnames = [
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
        
class HierarchicalSystem(System):
    
    def __init__(self, o1, o2):
        super().__init__(o1, o2)
        
        if not isinstance(o1, PrimaryOrbit):
            msg = ("o1 should be a PrimaryOrbit")
            raise AttributeError(msg)
        if not isinstance(o2, SatelliteOrbit):
            msg = ("o2 should be a SatelliteOrbit")
            raise AttributeError(msg)
            
        self.coordfunc = kep.coords_sat
        self.impfunc = kep.impacts_sat
        self.gimpfunc = kep.grad_impacts_sat  
        
        self.argnames = [
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
        