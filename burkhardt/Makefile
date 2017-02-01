
CXX = $(HOME)/bin/bin/g++ -std=c++17 -Wno-psabi

INC_DIR = ..
OBJ_DIR = obj

SRCS = \
  asa109.cpp \
  asa243.cpp \
  asa310.cpp \
  bernstein_polynomial.cpp \
  cdflib.cpp \
  chebyshev_polynomial.cpp \
  gegenbauer_polynomial.cpp \
  hermite_polynomial.cpp \
  jacobi_polynomial.cpp \
  kronrod.cpp \
  lobatto_polynomial.cpp \
  minmax.cpp \
  polpak.cpp \
  quadrule.cpp \
  test_values.cpp \
  timestamp.cpp \
  toms462.cpp \
  quadrule.cpp

OBJS = \
  $(OBJ_DIR)/asa109.o \
  $(OBJ_DIR)/asa243.o \
  $(OBJ_DIR)/asa310.o \
  $(OBJ_DIR)/bernstein_polynomial.o \
  $(OBJ_DIR)/cdflib.o \
  $(OBJ_DIR)/chebyshev_polynomial.o \
  $(OBJ_DIR)/gegenbauer_polynomial.o \
  $(OBJ_DIR)/hermite_polynomial.o \
  $(OBJ_DIR)/jacobi_polynomial.o \
  $(OBJ_DIR)/kronrod.o \
  $(OBJ_DIR)/lobatto_polynomial.o \
  $(OBJ_DIR)/minmax.o \
  $(OBJ_DIR)/polpak.o \
  $(OBJ_DIR)/timestamp.o \
  $(OBJ_DIR)/toms462.o \
  $(OBJ_DIR)/quadrule.o


libburkhardt.so: $(OBJS)
	$(CXX) -shared -o libburkhardt.so $(OBJS)
	cp libburkhardt.so libburkhardt.dll
	cp libburkhardt.so libburkhardt.dll ..

$(OBJ_DIR)/asa109.o: asa109.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asa109.o asa109.cpp

$(OBJ_DIR)/asa243.o: asa243.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asa243.o asa243.cpp

$(OBJ_DIR)/asa310.o: asa310.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asa310.o asa310.cpp

$(OBJ_DIR)/bernstein_polynomial.o: bernstein_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/bernstein_polynomial.o bernstein_polynomial.cpp

$(OBJ_DIR)/cdflib.o: cdflib.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cdflib.o cdflib.cpp

$(OBJ_DIR)/chebyshev_polynomial.o: chebyshev_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chebyshev_polynomial.o chebyshev_polynomial.cpp

$(OBJ_DIR)/gegenbauer_polynomial.o: gegenbauer_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gegenbauer_polynomial.o gegenbauer_polynomial.cpp

$(OBJ_DIR)/hermite_polynomial.o: hermite_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hermite_polynomial.o hermite_polynomial.cpp

$(OBJ_DIR)/jacobi_polynomial.o: jacobi_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jacobi_polynomial.o jacobi_polynomial.cpp

$(OBJ_DIR)/kronrod.o: kronrod.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/kronrod.o kronrod.cpp

$(OBJ_DIR)/lobatto_polynomial.o: lobatto_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/lobatto_polynomial.o lobatto_polynomial.cpp

$(OBJ_DIR)/minmax.o: minmax.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/minmax.o minmax.cpp

$(OBJ_DIR)/polpak.o: polpak.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polpak.o polpak.cpp

$(OBJ_DIR)/test_values.o: test_values.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/test_values.o test_values.cpp

$(OBJ_DIR)/timestamp.o: timestamp.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/timestamp.o timestamp.cpp

$(OBJ_DIR)/toms462.o: toms462.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/toms462.o toms462.cpp

$(OBJ_DIR)/quadrule.o: quadrule.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/quadrule.o quadrule.cpp

