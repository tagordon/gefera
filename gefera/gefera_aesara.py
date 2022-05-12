import aesara
from ops import ImpactsOp, FluxOp

__args__ = [
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

def validate_elements(argdict):
    
    if set(argdict.keys()) == set(__args__):
        
        idm = {
            v: i 
            for i, v in enumerate(__args__)
        }
        
        args = tuple(
            dict(
                sorted(
                    argdict.items(), 
                    key=lambda p: idm[p[0]]
                )
            ).values()
        )
        return args
        
    else:
        raise ValueError(
            "required parameters are: ", 
            __args__
        )
    
class BarycenterOrbit:
    
    """
    The orbit of the planet-moon barycenter around the star. 
    
    Args:
        ab (scalar): Semimajor axis
        tb (scalar): Time of periastron passage
        eb (scalar): Eccentricity
        pb (scalar): Period
        wb (scalar): Argument of periastron (in radians)
        ib (scalar): Inclination (in radians)
    """
    
    def __init__(self, ab, tb, eb, pb, wb, ib):
        
        self.pnames = [
            'ab', 
            'tb', 
            'eb', 
            'pb', 
            'wb', 
            'ib'
        ]
        self.pdict = dict(
            zip(
                self.pnames,
                (ab, tb, eb, pb, wb, ib)
            )
        )
    
class MoonOrbit:
    
    """
    The orbit of the moon around the planet.
    
    Args:
        am (scalar): Semimajor axis
        tm (scalar): Time of periastron passage
        em (scalar): Eccentricity
        pm (scalar): Period
        om (scalar): Longitude of ascending node (in radians)
        wm (scalar): Argument of periastron (in radians)
        im (scalar): Inclination (in radians)
        mm (scalar): Moon/planet mass ratio
    """
    
    def __init__(self, am, tm, em, pm, om, wm, im, mm):
        
        self.pnames = [
            'am', 
            'tm', 
            'em', 
            'pm', 
            'om', 
            'wm', 
            'im', 
            'mm'
        ]
        self.pdict = dict(
            zip(
                self.pnames, 
                (am, tm, em, pm, om, wm, im, mm)
            )
        )
        
class System:
    
    """
    Class representing the three-body star-planet-moon system. 
    
    
    .. note:: The dynamics of the system are approximated by assuming 
    no gravitational interaction between the moon and the star. 
    In other words, the planet/moon system is treated as a single 
    body for purposes of computing the orbit of their mutual 
    barycenter around the star, and then the orbit of the moon 
    around the planet is computed without considering star-moon 
    interactions. This apprximation should be valid in the limit 
    that the moon's orbit is well within the planet's Hill sphere. 
    
    Args:
        bo (BarycenterOrbit): The orbit of the planet-moon 
                barycenter with respect to the star.
        mo (MoonOrbit): The orbit of the moon with respect to the 
                planet. 
    """
    
    def __init__(self, bo, mo):
        
        self.bo = bo
        self.mo = mo
        
    def lightcurve(self, t, u1, u2, rp, rm):
        
        """
        Get the lightcurve resulting from a transit of the moon/planet system.
        
        Args: 
            t (tensor): Times at which the flux should be computed
            rp (scalar): Radius of the planet
            rm (scalar): Radius of the moon
            u1 (scalar): The first limb-darkening parameter
            u2 (scalar): The second limb-darkening parameter
        """
        
        impacts = ImpactsOp(t)
        flux = FluxOp()
        
        impact_params = validate_elements(
            {
                **self.bo.pdict, 
                **self.mo.pdict
            }
        )
        
        out = impacts(*impact_params)
        bp, bpm, theta = out[0], out[1], out[2]
        return flux(
            u1, 
            u2, 
            rp, 
            rm, 
            bp, bpm, theta
        )[0]