# This Makefile builds a test program in Fortran for lerchphi.c 
# in standard double precision.
#
# Program lerchphi.c is copyright by
#
# Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
# Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.
#
# Version 1.00 (May 1, 2002)
#
#
# Please set the following variables:
#
# Set FC to whatever Fortran compiler is on your system (e.g., f77 or g77)
# Set CC to whatever C compiler is on your system (e.g., cc or gcc)
# Set LIB to the standard system math library (e.g., lm)
# Set INCPATH to where your headers are
# Set LIBPATH to where your libraries are
# Set PREP to ADD_UNDERSCORE for Sun Fortran compiler

FC = f77
CC = cc
LIB = lm
INCPATH = ./
LIBPATH = ./
PREP = -DADD_UNDERSCORE

testf77: testf77.o lerchphi.f77.o
	$(FC) -L$(LIBPATH) -o testf77 testf77.o lerchphi.f77.o -$(LIB)

testf77.o: testf77.f
	$(FC) -I$(INCPATH) -c testf77.f

lerchphi.f77.o: ../Source/lerchphi.c
	$(CC) $(PREP) -I$(INCPATH) -o lerchphi.f77.o -c ../Source/lerchphi.c

clean:
	rm testf77 testf77.o lerchphi.f77.o
