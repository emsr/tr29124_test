
//SUFFIX = _tr29124
SUFFIX = _specfun

CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
CXX_SRC_DIR = $(HOME)/gcc$(SUFFIX)

GCC = $(CXX_INST_DIR)/bin/gcc -g -Wall -Wextra
CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi -I..
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/7.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64

OBJ_DIR = obj

BINS = hermite_test \
       laguerre_test \
       legendre_test \
       gegenbauer_test \
       jacobi_test \
       chebyshev_t_test \
       chebyshev_u_test \
       chebyshev_v_test \
       chebyshev_w_test \
       radpoly_test \
       zernike_test


all: $(BINS)


hermite_test: $(OBJ_DIR)/hermite_test.o
	$(CXX17) -o hermite_test $(OBJ_DIR)/hermite_test.o -lquadmath

laguerre_test: $(OBJ_DIR)/laguerre_test.o
	$(CXX17) -o laguerre_test $(OBJ_DIR)/laguerre_test.o -lquadmath

legendre_test: $(OBJ_DIR)/legendre_test.o
	$(CXX17) -o legendre_test $(OBJ_DIR)/legendre_test.o -lquadmath

gegenbauer_test: $(OBJ_DIR)/gegenbauer_test.o
	$(CXX17) -o gegenbauer_test $(OBJ_DIR)/gegenbauer_test.o -lquadmath

jacobi_test: $(OBJ_DIR)/jacobi_test.o
	$(CXX17) -o jacobi_test $(OBJ_DIR)/jacobi_test.o -lquadmath

chebyshev_t_test: $(OBJ_DIR)/chebyshev_t_test.o
	$(CXX17) -o chebyshev_t_test $(OBJ_DIR)/chebyshev_t_test.o -lquadmath

chebyshev_u_test: $(OBJ_DIR)/chebyshev_u_test.o
	$(CXX17) -o chebyshev_u_test $(OBJ_DIR)/chebyshev_u_test.o -lquadmath

chebyshev_v_test: $(OBJ_DIR)/chebyshev_v_test.o
	$(CXX17) -o chebyshev_v_test $(OBJ_DIR)/chebyshev_v_test.o -lquadmath

chebyshev_w_test: $(OBJ_DIR)/chebyshev_w_test.o
	$(CXX17) -o chebyshev_w_test $(OBJ_DIR)/chebyshev_w_test.o -lquadmath

radpoly_test: $(OBJ_DIR)/radpoly_test.o
	$(CXX17) -o radpoly_test $(OBJ_DIR)/radpoly_test.o -lquadmath

zernicke_test: $(OBJ_DIR)/zernicke_test.o
	$(CXX17) -o zernicke_test $(OBJ_DIR)/zernicke_test.o -lquadmath


$(OBJ_DIR)/hermite_test.o: *.h hermite_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/hermite_test.o hermite_test.cpp

$(OBJ_DIR)/laguerre_test.o: *.h laguerre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/laguerre_test.o laguerre_test.cpp

$(OBJ_DIR)/legendre_test.o: *.h legendre_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/legendre_test.o legendre_test.cpp

$(OBJ_DIR)/gegenbauer_test.o: *.h gegenbauer_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/gegenbauer_test.o gegenbauer_test.cpp

$(OBJ_DIR)/jacobi_test.o: *.h jacobi_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/jacobi_test.o jacobi_test.cpp

$(OBJ_DIR)/chebyshev_t_test.o: *.h chebyshev_t_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_t_test.o chebyshev_t_test.cpp

$(OBJ_DIR)/chebyshev_u_test.o: *.h chebyshev_u_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_u_test.o chebyshev_u_test.cpp

$(OBJ_DIR)/chebyshev_v_test.o: *.h chebyshev_v_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_v_test.o chebyshev_v_test.cpp

$(OBJ_DIR)/chebyshev_w_test.o: *.h chebyshev_w_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/chebyshev_w_test.o chebyshev_w_test.cpp

$(OBJ_DIR)/radpoly_test.o: *.h radpoly_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/radpoly_test.o radpoly_test.cpp

$(OBJ_DIR)/zernicke_test.o: *.h zernicke_test.cpp
	$(CXX17) -c -o $(OBJ_DIR)/zernicke_test.o zernicke_test.cpp

