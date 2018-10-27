# This Makefile builds a test program for lerchphi.c 
# in standard double precision.
#
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
# Set CC to whatever C compiler is on your system (e.g., cc or gcc)
# Set LIB to the standard system math library (e.g., -lm)
# Set INCPATH to where your headers are
# Set LIBPATH to where your libraries are

CC = cc
LIB = lm
INCPATH = ./
LIBPATH = ./

test : test.o lerchphi.o
	$(CC) -L$(LIBPATH) -o test test.o lerchphi.o -$(LIB)

test.o : test.c
	$(CC) -I$(INCPATH) -g -c test.c
lerchphi.o : ../Source/lerchphi.c
	$(CC) -I$(INCPATH) -g -c ../Source/lerchphi.c

clean : 
	rm test test.o lerchphi.o
