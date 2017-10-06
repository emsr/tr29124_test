all: test_rpoly test_bairstow test_broyden test_broyden_cmplx test_newton test_newton_cmplx

test_rpoly: rpoly.cpp
	g++ -g -o test_rpoly rpoly.cpp

test_bairstow: bairstow.c
	gcc -g -o test_bairstow bairstow.c

test_broyden: broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp
	g++ -g -o test_broyden broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp

test_broyden_cmplx: broyden_cmplx/cbroyden.cpp broyden_cmplx/cbtst.cpp broyden_cmplx/cged.cpp
	g++ -g -o test_broyden_cmplx broyden_cmplx/cbroyden.cpp broyden_cmplx/cbtst.cpp broyden_cmplx/cged.cpp

test_newton: newton/newt.h newton/ged.cpp newton/nlnewt.cpp
	g++ -g -o test_newton newton/ged.cpp newton/nlnewt.cpp

test_newton_cmplx: newton_cmplx/newt.h newton_cmplx/cged.cpp newton_cmplx/cnlnewt.cpp
	g++ -g -o test_newton_cmplx newton_cmplx/cged.cpp newton_cmplx/cnlnewt.cpp

