# Makefile for the double delta function integral with a weight in the BZ

FFTW3_HOME="${HOME}/usr"
# gfortran
FC = gfortran
#LIB=-I${FFTW3_HOME}/include -L${FFTW3_HOME}/lib -lfftw3 [-static]
LIB= -I/usr/include -L/usr/lib -lfftw3
FFLAGS=-g -Wall -ffree-form -O2 -fopenmp 

test: declarations.o fftw.o compute.o readio.o dos.o
	$(FC) $(FFLAGS) $(LIB) -o test declarations.o fftw.o compute.o readio.o dos.o



declarations.mod: declarations.o declarations.f90
	$(FC) $(FFLAGS) -c declarations.f90
declarations.o: declarations.f90
	$(FC) $(FFLAGS) -c declarations.f90
fftw.mod : fftw.o fftw.f90
	$(FC) $(FFLAGS) $(LIB) -c fftw.f90
fftw.o : fftw.f90
	$(FC) $(FFLAGS) $(LIB) -c fftw.f90
readio.mod : readio.o readio.f90
	$(FC) $(FFLAGS) -c readio.f90
readio.o : declarations.mod readio.f90
	$(FC) $(FFLAGS) -c readio.f90
compute.mod: compute.o compute.f90
	$(FC) $(FFLAGS) -c compute.f90
compute.o: declarations.mod fftw.mod compute.f90
	$(FC) $(FFLAGS) -c compute.f90
dos.o: declarations.mod fftw.mod readio.mod compute.mod dos.f90
	$(FC) $(FFLAGS) -c dos.f90

clean:
	rm *.mod *.o
