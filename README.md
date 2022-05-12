# gefera

## installation

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

Now cd into the build directory and issue the command `meson compile`. This will produce the shared fortran libraries that gefera is built on. Once the compilation is finished return to the main directory and use pip to install the gefera. 

```
cd ../; pip install .
```
 
## basic usage

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
