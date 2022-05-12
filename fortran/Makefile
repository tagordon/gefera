all: ellip.o kepler phot.o
	gfortran -shared -Ofast ellip.o phot.o phot_nograd.o -o wrapper.so
	gfortran -shared -Ofast kep.o hrch.o -o hrch.so
	gfortran -shared -Ofast kep.o conf.o -o conf.so

ellip.o: 
	gfortran -g -shared -Ofast -unroll=12 -c ellip.f90

kepler: kep.o
	gfortran -g -shared -Ofast -unroll=12 -fmax-stack-var-size=100000 -c conf.f90
	gfortran -g -shared -Ofast -unroll=12 -fmax-stack-var-size=100000 -c hrch.f90

kep.o:
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 kep.f90
phot.o: 
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 phot.f90
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 phot_nograd.f90
clean: 
	rm *.o
	rm *.so
	rm *.mod
