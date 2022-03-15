
#######################################################################
#                                                                     #
#             Makefile for ORTHPOL                                    #
#                                                                     #
#  Usage:                                                             #
#                                                                     #
#  Type  make                  to compile and run all test programs   #
#  Type  make testN.out.local  to compile and run the test program    #
#                              number N                               #
#  Type  make clean            to remove all but the generated output #
#                                                                     #
#######################################################################

#----------------------------------------------------------------------
# This section specifies the name of the linker and compilation options
# for various systems. All but one are commented out.
#----------------------------------------------------------------------

# ... For Linux
LD = gfortran
#FFLAGS =  #Do not use -O2; it may introduce errors

# ... For the Sun
#LD = f77
#FFLAGS =  #Do not use -O2; it may introduce errors

# ... For the Cray
#LD = cf77
#FFLAGS =

# ... For the Convex
#LD = fc
#FFLAGS = -O2

# ... For the RS/6000
#LD = xlf
#FFLAGS = -O

# ... For the HP 9000 series
#LD = f77
#FFLAGS = -O -K -w

# ... For Silicon Graphics
#LD = f77
#FFLAGS = -O2

#-------------------------------------------------------------------
# This section specifies which object files are needed for each test
# program.
#-------------------------------------------------------------------

TEST1 = test1.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o
TEST2 = test2.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o
TEST3 = test3.o r1mach.o sti.o lancz.o
TEST4 = test4.o r1mach.o mcdis.o qgp.o sti.o lancz.o recur.o gauss.o
TEST5 = test5.o r1mach.o d1mach.o recur.o drecur.o mcdis.o qgp.o sti.o lancz.o gauss.o
TEST6 = test6.o r1mach.o d1mach.o mcdis.o dmcdis.o qgp.o dqgp.o sti.o dsti.o lancz.o dlancz.o recur.o drecur.o gauss.o dgauss.o
TEST7 = test7.o r1mach.o d1mach.o mcdis.o dmcdis.o qgp.o dqgp.o sti.o dsti.o lancz.o dlancz.o
TEST8 = test8.o r1mach.o recur.o mccheb.o qgp.o cheb.o
TEST9 = test9.o r1mach.o d1mach.o recur.o drecur.o cheb.o dcheb.o chri.o dchri.o
TEST10= test10.o r1mach.o d1mach.o recur.o drecur.o gauss.o dgauss.o chri.o dchri.o
TEST11= test11.o r1mach.o d1mach.o recur.o drecur.o knum.o nu0jac.o dknum.o gchri.o dgchri.o cheb.o dcheb.o chri.o dchri.o

#--------------------------------------------------------------
# This is the default target. It will cause all the tests to be
# generated and run, saving their outputs in files.
#--------------------------------------------------------------

all :   test1.out.local test2.out.local test3.out.local test4.out.local test5.out.local test6.out.local test7.out.local test8.out.local test9.out.local test10.out.local test11.out.local

#----------------------------------------------------------------------
# This section shows how to run each test program. The output for testN
# is put in testN.out.local.
#----------------------------------------------------------------------

test1.out.local  : test1  ; ./test1  > $@
test2.out.local  : test2  ; ./test2  > $@
test3.out.local  : test3  ; ./test3  > $@
test4.out.local  : test4  ; ./test4  > $@
test5.out.local  : test5  ; ./test5  > $@
test6.out.local  : test6  ; ./test6  > $@
test7.out.local  : test7  ; ./test7  > $@
test8.out.local  : test8  ; ./test8  > $@
test9.out.local  : test9  ; ./test9  > $@
test10.out.local : test10 ; ./test10 > $@
test11.out.local : test11 ; ./test11 > $@

#--------------------------------------------------
# This section shows how to link each test program.
#--------------------------------------------------

test1  : $(TEST1)  ; $(LD) -o $@ $(TEST1)
test2  : $(TEST2)  ; $(LD) -o $@ $(TEST2)
test3  : $(TEST3)  ; $(LD) -o $@ $(TEST3)
test4  : $(TEST4)  ; $(LD) -o $@ $(TEST4)
test5  : $(TEST5)  ; $(LD) -o $@ $(TEST5)
test6  : $(TEST6)  ; $(LD) -o $@ $(TEST6)
test7  : $(TEST7)  ; $(LD) -o $@ $(TEST7)
test8  : $(TEST8)  ; $(LD) -o $@ $(TEST8)
test9  : $(TEST9)  ; $(LD) -o $@ $(TEST9)
test10 : $(TEST10) ; $(LD) -o $@ $(TEST10)
test11 : $(TEST11) ; $(LD) -o $@ $(TEST11)

#----------------------------------------------
# This section shows how to clean up afterward.
#----------------------------------------------

clean :
	rm -f core test1 test2 test3 test4 test5 test6 test7 test8 test9 test10 test11 *.o

