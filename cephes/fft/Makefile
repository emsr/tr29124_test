#				Makefile
# 	to build and test the Numerical Math Library (libla.a)
#
# To build the library, you may want to edit the list of library modules
# below (MODULES). Say, if you don't need/want FFT, remove all the names
# that contain fft from MODULES=. Then point ./c++ to the right
# script: g++.2.7.2 or g++.2.8.1 (see below if using a non-GNU C++ compiler).
# Then make the modules you chose:
#		make lib
# A simple 'make' would suffice, too.
#
# To verify the library, do
#	make check-all
# or, specifically,
#	make vmatrix		(checks elementary matrix operations)
#	make vmatrix1		(checks less trivial matrix operations)
#	make vmatrix2		(checks advanced matrix operations: multipl)
#	make vlastreams		(checks treams over matrices and vectors)
#	make vvector		(checks elementary vector operations)
#	
#	make vzeroin		(verifies Brent's 1D root finder)
#	make vfminbr		(verifies Brent's 1D minimizer)
#
#	make vali		(check Aitken-Langrange interpolators)
#	make vhjmin		(check a multi-variate optimizer)
#
#	make vsvd		(check SVD decompositions)
#	make vslesing		(test solving generalized Ax=b)
#
#	make vfft		(check and clock FFT)
#
#	make sample_ult
#	make sample_adv
# Note, this Makefile was built (and works under) GNU make 3.71
#
# Please check RANLIB below and adjust it for your system if necessary
# (it was made for a BSD-like system)
# $Id: Makefile,v 4.3 1998/12/25 23:05:42 oleg Exp oleg $

# Make sure ./c++ points to the right g++ script.
# If using a non-GNU compiler, setting
# CC=cc -c  and CCL=cc is a good starting point
CC=./c++
CCL=g++ -pipe
.SUFFIXES: .cc
MODULES=myenv.cc matrix1.cc matrix2.cc matrix_sub.cc \
	vector.cc determinant.cc matrix_inv.cc \
	zeroin.cc fminbr.cc ali.cc hjmin.cc svd.cc \
	fft_init.cc fft_input.cc fft_output.cc
LIBRARY=libla.a
CLDEF=			# Other flags for the linker, if needed
LIBS=-lm		# Other libraries to use, if needed \
			# When compiling FFT with SunWPro, add -lcomplex
#RANLIB = (ar d $(LIBRARY) __.SYMDEF || true); ranlib $(LIBRARY) # BSD-like
RANLIB=/bin/true		# SYSV-based
#	Rules, new style

%.o : %.cc
	$(CC) $*.cc

$(LIBRARY)(%.o) : %.cc
	$(CC) $*.cc
	ar rcv $(LIBRARY) $*.o
	rm -f $*.o

% : %.o $(LIBRARY)
	$(CCL) $< $(CLDEF) $(LIBRARY) $(LIBS) -o $@
	./$@

% :: %.cc
	$(CC) $*.cc
	$(CCL) $*.o $(CLDEF) $(LIBRARY) $(LIBS) -o $@
	./$@

#	Rules, old style
#.o:	$*.o $(LIBRARY)
#	$(CCL) $*.o $(CLDEF) $(LIBRARY) $(LIBS) -o $*
#	./$*
#.cc: 	$*.cc $(LIBRARY)
#	$(CC) $*.cc
#	$(CCL) $*.o $(CLDEF) $(LIBRARY) $(LIBS) -o $*
#	./$*
#.cc.o:
#	$(CC) $*.cc
#

# Primary goal

# Library

lib:	$(LIBRARY)
.PHONY:		lib
.PRECIOUS:	$(LIBRARY)

# 			Compile the source files that have been changed 

$(LIBRARY):    $(LIBRARY)($(MODULES:.cc=.o))
	$(RANLIB)

#$(LIBRARY)::    $(MODULES)
#	$(CC) $?
#	listobj=`echo $? | sed s/.cc/.o/g` ; \
#	ar rv $(LIBRARY) $$listobj &&	\
#	rm $$listobj
#	$(RANLIB)

# Verification routines
check-all:	vmatrix vvector vmatrix1 vmatrix2 vlastreams \
		vali vhjmin vfminbr vzeroin \
		vsvd vslesing vfft

clean:
	rm -f *.o core vmatrix vvector vmatrix1 vmatrix2 vlastreams \
		vali vhjmin \
		vfminbr vzeroin sample_ult sample_adv \
		vsvd vslesing vfft

dist-clean:	clean
	rm -f $(LIBRARY)


# Specific dependent goals

# Optimization causes g++'s internal compiler error...
#$(LIBRARY)(matrix1.o): matrix1.cc
#	$(CC) -O0 matrix1.cc
#	ar rcv $(LIBRARY) matrix1.o
#	rm -f matrix1.o
#
#vmatrix.o : vmatrix.cc
#	$(CC) -O0 vmatrix.cc

# Dependence rules

matrix1.cc:	LinAlg.h

#$(LIBRARY)::	LinAlg.h
#	$(MAKE) -W matrix1.cc lib

