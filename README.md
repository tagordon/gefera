![tests](https://github.com/tagordon/gefera/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/gefera/badge/?version=latest)](https://gefera.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/491390764.svg)](https://zenodo.org/badge/latestdoi/491390764)

# gefera

[Installation](#installation)

[Basic Usage](#basic-usage)

[The Likelihood](#the-likelihood)

## Installation

To install gefera run 

```
pip install gefera
```

Or clone the github repository and install from the source code using pip: 

```
git clone https://www.github.com/tagordon/gefera
cd gefera; pip install .
``` 
 
## Basic Usage

Gefera has two components: the dynamical model which is used to compute the positions of two bodies at a given time or array of times, and the photometric model which computes the observed flux given the positions of the bodies. Gefera implements two dynamical systems. The first is a `HierarchicalSystem`, which describes an exomoon/exoplanet pair. To instantiate this system let's start by defining the dynamical parameters:

For the planet:

```python
ap = 215.0			# semimajor axis in stellar radii
tp = -91.25			# starting epoch in days
ep = 0.0			# eccentricity 
pp = 365			# period in days
wp = 0.1 * np.pi / 180		# longitude of periastron in degrees
ip = 89.8 * np.pi / 180		# inclination in degrees
```

and for the moon:    

```python
am = 2			# semimajor axis of the moon's orbit around the planet in stellar radii
tm = -4.2			# starting epoch in days
em = 0.0			# eccentricity
pm = 8				# period in days
om = 90 * np.pi / 180		# longitude of the ascending node in degrees
wm = -90 * np.pi / 180		# longitude of periastron in degrees
im = 90.0 * np.pi / 180		# inclination in degrees
mm = 0.01			# mass of the moon in units of the mass of the planet
```

We can then instantiate the system as follows: 

```python
import gefera as gf

po = gf.orbits.PrimaryOrbit(ap, tp, ep, pp, wp, ip)
mo = gf.orbits.SatelliteOrbit(am, tm, em, pm, om, wm, im, mm)
sys = gf.systems.HierarchicalSystem(po, mo)
```

To compute the light curve we also need the radii of both bodies and the quadratic limb darkening parameters (right now gefera only supports quadratic limb darkening)

```python
rp = 0.1
rm = 0.05
u1 = 0.5
u2 = 0.3
```

We can then compute the light curve:

```python
t = np.linspace(-0.6, 0.3, 10000)
# The out-of-transit flux is set to zero by default, so we add 1 to get the normalized flux.
flux = sys.lightcurve(t, u1, u2, rp, rm) + 1
```
When plotted, the lightcurve for this system looks like this:
![exomoon lightcurve](https://github.com/tagordon/gefera/blob/master/notebooks/readme_plot_light.png#gh-light-mode-only)
![exomoon lightcurve](https://github.com/tagordon/gefera/blob/master/notebooks/readme_plot_dark.png#gh-dark-mode-only)

If the gradient is required we use

```python
flux, grad = sys.lightcurve(t, u1, u2, rp, rm, grad=True)
```

This will return the gradient of the light curve with respect to each of the input parameters as a dictionary keyed by the name of the parameter. The keys are `a1, t1, e1... a2, t2, e2...` where the `1` parameters are for the larger body and `2` refers to the smaller body. 

## The Likelihood 

In order to conduct inference using the gefera model we need to compute the likelihood with respect to some observations. If `t` is an array containing the times of the observations and `y` is an array containing flux measurements (normalized to an out-of-transit flux of 1) then we can compute the likelihood for the model as follows:

```python
ll = sys.loglike(y - 1, t, u1, u2, r1, r2, sigma)
``` 

where `sigma` is the uncertainty, assumed to be the same for each measurement. If we want to make use of the gradient in our inference then we can call:

```python
ll, dll = sys.loglike(y - 1, t, u1, u2, r1, r2, sigma, grad=True)
```

The gradient is returned as an array of arrays containing the gradients with respect to each parameter in the order `a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, om2, w2, i2, mm2, r1, r2, u1, u2`. If `func` is a callable function that returns `-ll, -dll` and `init` is an array of initial parameter guesses in the same order as the gradients are returned then we can use `scipy.optimize.minimize` to optimize the model as follows:

```python
res = minimize(func, init, jac=True, method='TNC')
```
