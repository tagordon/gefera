all: ellip.o kep.o phot.o
	gfortran -shared -Ofast ellip.o phot.o phot_nograd.o kep.o -o wrapper.so

ellip.o: 
	gfortran -g -shared -Ofast -unroll=12 -c ellip.f90

kep.o:
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 kep.f90

phot.o: 
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 phot.f90
	gfortran -g -shared -Ofast -c -unroll=5 -fmax-stack-var-size=100000 phot_nograd.f90
clean: 
	rm *.o
	rm *.so
	rm *.mod
