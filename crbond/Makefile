all: test_bairstow test_broyden

test_bairstow: bairstow.c
	gcc -g -o test_bairstow bairstow.c

test_broyden: broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp
	g++ -g -o test_broyden broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp

#test_bairstow: test_bairstow.cpp
#	g++ -g -o test_bairstow test_bairstow.cpp

