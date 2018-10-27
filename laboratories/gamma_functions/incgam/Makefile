
FC = $(HOME)/bin/bin/gfortran

testincgam: cons.f90 incgam.f90 testincgam.f90
	$(FC) -c cons.f90
	$(FC) -c incgam.f90
	$(FC) -c testincgam.f90
	$(FC) cons.o incgam.o testincgam.o -o testincgam

test: testincgam
	LD_LIBRARY_PATH=$(HOME)/bin/lib64:$(LD_LIBRARY_PATH) ./testincgam

