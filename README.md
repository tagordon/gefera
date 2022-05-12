# gefera

[Installation](#installation)

[Basic Usage](#installation)

[The Likelihood](#the-likelihood)

## Installation

Gefera uses the Meson build system to compile the fortran backend. If you haven't used Meson before, you'll need to start by installing it via pip: 

```
pip install meson
```

To install gefera begin by cloning the github repository. 

```
git clone https://www.github.com/tagordon/gefera
```

Once the repository has downloaded, cd into the main directory and run `meson setup builddir` to make the build directory.

```
cd gefera; meson setup builddir
``` 

Now cd into the build directory and issue the command `meson compile`:

```
cd builddir; meson compile
```

This will produce the shared fortran libraries that gefera is built on. Once the compilation is finished return to the main directory and use pip to install the gefera. 

```
cd ../; pip install .
```
 
## Basic Usage

Gefera has two components: the dynamical model which is used to compute the positions of two bodies at a given time or array of times, and the photometric model which computes the observed flux given the positions of the bodies. Gefera implements two dynamical systems. The first is a `HierarchicalSystem`, which describes an exomoon/exoplanet pair. To instantiate this system let's start by defining the dynamical parameters:

For the planet:

```python
ap = 1.0			# semimajor axis in au
tp = -91.25			# starting epoch in days
ep = 0.0			# eccentricity 
pp = 365			# period in days
wp = 0.1 * np.pi / 180		# longitude of periastron in degrees
ip = 89.8 * np.pi / 180		# inclination in degrees
```

and for the moon:    

```python
am = 0.02			# semimajor axis of the moon's orbit around the planet in au
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
po = gf.PrimaryOrbit(ap, tp, ep, pp, wp, ip)
mo = gf.SatelliteOrbit(am, tm, em, pm, om, wm, im, mm)
sys = gf.HierarchicalSystem(po, mo)
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
flux = sys.lightcurve(t, u1, u2, rp, rm)
```

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
ll, dll = sys.loglike(y - 1, t, u1, u2, r1, r2, sigma, grad=True, sign=-1)
```

Here `sign=-1` indicates that we want the negative of the log-likelihood and its gradient. If we want the positive log-likelihood and gradient we can omit this argument. The gradient is returned as an array of arrays containing the gradients with respect to each parameter in the order `a1, t1, e1, p1, w1, i1, a2, t2, e2, p2, om2, w2, i2, mm2, r1, r2, u1, u2`. If `func` is a callable function that returns `ll, dll` and `init` is an array of initial parameter guesses in the same order as the gradients are returned then we can use `scipy.optimize.minimize` to optimize the model as follows:

```python
res = minimize(func, init, jac=True, method='TNC')
```
