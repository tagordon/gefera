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
 
## usage
