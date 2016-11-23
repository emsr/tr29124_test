
CXX = $(HOME)/bin/bin/g++ -std=c++17 -Wno-psabi

INC_DIR = ..
OBJ_DIR = obj

SRCS = \
  asa109.cpp \
  asa243.cpp \
  asa310.cpp \
  bernstein_polynomial.cpp \
  cdflib.cpp \
  hermite_polynomial.cpp \
  lobatto_polynomial.cpp \
  minmax.cpp \
  polpak.cpp \
  test_values.cpp \
  timestamp.cpp \
  toms462.cpp

OBJS = \
  $(OBJ_DIR)/asa109.o \
  $(OBJ_DIR)/asa243.o \
  $(OBJ_DIR)/asa310.o \
  $(OBJ_DIR)/bernstein_polynomial.o \
  $(OBJ_DIR)/cdflib.o \
  $(OBJ_DIR)/hermite_polynomial.o \
  $(OBJ_DIR)/lobatto_polynomial.o \
  $(OBJ_DIR)/minmax.o \
  $(OBJ_DIR)/polpak.o \
  $(OBJ_DIR)/test_values.o \
  $(OBJ_DIR)/timestamp.o \
  $(OBJ_DIR)/toms462.o

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

$(OBJ_DIR)/hermite_polynomial.o: hermite_polynomial.cpp
	$(CXX) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hermite_polynomial.o hermite_polynomial.cpp

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

