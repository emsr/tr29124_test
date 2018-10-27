# This Makefile builds a test program for lerchphimp.c 
# in quad-double precision.  It needs quad-double library
# by David M. Bailey.
#
#
# Programs lerchphi.c and lerchphimp.cc are copyright by
#
# Sergej V. Aksenov (http://www.geocities.com/saksenov) and 
# Ulrich D. Jentschura (jentschura@physik.tu-dresden.de), 2002.
#
# Version 1.00 (May 1, 2002)
#
#
# Please set the following variables:
#
# Set CC to whatever C++ compiler is on your system (e.g., c++ or g++)
# Set LIB to the quad double library (e.g., -lqd)
# Set INCPATH to where your headers are
# Set LIBPATH to where your libraries are
# Set PREP to any options for the compiler.  In particular, set 
# the same options used to compile libqd library.

CC = c++
LIB = lqd
INCPATH = ../../qd/include
LIBPATH = ../../qd
PREP = -DACCURATE # -Dx86 -Wall -O

testmp : testmp.o lerchphimp.o
	$(CC) $(PREP) -L$(LIBPATH) -o testmp testmp.o lerchphimp.o -$(LIB)

testmp.o : testmp.cc
	$(CC) $(PREP) -I$(INCPATH) -g -c testmp.cc
lerchphimp.o : ../Source/lerchphimp.cc
	$(CC) $(PREP) -I$(INCPATH) -g -c ../Source/lerchphimp.cc

clean : 
	rm testmp testmp.o lerchphimp.o
