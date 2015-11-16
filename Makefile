
#HOME = /home/ed
# -Wconversion
CXX = $(HOME)/bin_specfun/bin/g++ -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__
CXX_INC_DIR = $(HOME)/bin_specfun/include/c++/6.0.0/bits
CXX_LIB_DIR = $(HOME)/bin_specfun/lib64

BINS = test_special_function \
       diff_special_function \
       testcase \
       test_airy \
       test_csint \
       test_fresnel \
       test_hermite \
       test_limits \
       test_cmath \
       test_bessel \
       test_nric_bessel \
       test_legendre

all: test_special_function \
     diff_special_function \
     testcase \
     test_airy \
     test_csint \
     test_fresnel \
     test_hermite \
     test_limits \
     test_cmath \
     test_bessel \
     test_nric_bessel \
     test_legendre

testcases: testcase
	LD_LIBRARY_PATH=/home/ed/bin_specfun/lib64:$LD_LIBRARY_PATH ./testcase

test:
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_airy > test_airy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_csint > test_csint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_fresnel > test_fresnel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hermite > test_hermite.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_limits > test_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_cmath > test_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bessel > test_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_nric_bessel > test_nric_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./diff_special_function > diff_special_function.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_special_function > test_special_function.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_legendre > test_legendre.txt

check: \
       check_airy_ai \
       check_airy_bi \
       check_assoc_laguerre \
       check_assoc_legendre \
       check_beta \
       check_comp_ellint_1 \
       check_comp_ellint_2 \
       check_comp_ellint_3 \
       check_conf_hyperg \
       check_cyl_bessel_i \
       check_cyl_bessel_j \
       check_cyl_bessel_j_asymp \
       check_cyl_bessel_k \
       check_cyl_neumann \
       check_cyl_neumann_asymp \
       check_ellint_1 \
       check_ellint_2 \
       check_ellint_3 \
       check_expint_neg \
       check_expint_pos \
       check_hyperg \
       check_laguerre \
       check_legendre \
       check_riemann_zeta_neg \
       check_riemann_zeta_pos \
       check_sph_bessel \
       check_sph_legendre \
       check_sph_neumann



test_special_function: test_special_function.cpp gsl_wrap.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o test_special_function test_special_function.cpp gsl_wrap.cpp -lgsl -lgslcblas

diff_special_function: diff_special_function.cpp gsl_wrap.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o diff_special_function diff_special_function.cpp gsl_wrap.cpp -lgsl -lgslcblas

testcase: testcase.cpp testcase.tcc gsl_wrap.cpp $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o testcase testcase.cpp gsl_wrap.cpp -lgslcblas -lgsl

test_limits: test_limits.cpp
	$(CXX) -o test_limits test_limits.cpp

test_cmath: test_cmath.cpp
	$(CXX) -o test_cmath test_cmath.cpp

test_airy: test_airy.cpp airy.tcc
	$(CXX) -o test_airy test_airy.cpp

test_csint: test_csint.cpp csint.tcc
	$(CXX) -o test_csint test_csint.cpp

test_fresnel: test_fresnel.cpp fresnel.tcc
	$(CXX) -o test_fresnel test_fresnel.cpp

test_hermite: test_hermite.cpp new_hermite.tcc
	$(CXX) -o test_hermite test_hermite.cpp

test_bessel: test_bessel.cpp new_bessel.tcc
	$(CXX) -o test_bessel test_bessel.cpp

test_nric_bessel: test_nric_bessel.cpp nric_bessel.tcc
	$(CXX) -o test_nric_bessel test_nric_bessel.cpp

test_legendre: test_legendre.cpp legendre.tcc
	$(CXX) -o test_legendre test_legendre.cpp


check_airy_ai: check_airy_ai.cc
	$(CXX) -D__TEST_DEBUG -o check_airy_ai check_airy_ai.cc

check_airy_bi: check_airy_bi.cc
	$(CXX) -D__TEST_DEBUG -o check_airy_bi check_airy_bi.cc

check_assoc_laguerre: check_assoc_laguerre.cc
	$(CXX) -D__TEST_DEBUG -o check_assoc_laguerre check_assoc_laguerre.cc

check_assoc_legendre: check_assoc_legendre.cc
	$(CXX) -D__TEST_DEBUG -o check_assoc_legendre check_assoc_legendre.cc

check_beta: check_beta.cc
	$(CXX) -D__TEST_DEBUG -o check_beta check_beta.cc

check_comp_ellint_1: check_comp_ellint_1.cc
	$(CXX) -D__TEST_DEBUG -o check_comp_ellint_1 check_comp_ellint_1.cc

check_comp_ellint_2: check_comp_ellint_2.cc
	$(CXX) -D__TEST_DEBUG -o check_comp_ellint_2 check_comp_ellint_2.cc

check_comp_ellint_3: check_comp_ellint_3.cc
	$(CXX) -D__TEST_DEBUG -o check_comp_ellint_3 check_comp_ellint_3.cc

check_conf_hyperg: check_conf_hyperg.cc
	$(CXX) -D__TEST_DEBUG -o check_conf_hyperg check_conf_hyperg.cc

check_cyl_bessel_i: check_cyl_bessel_i.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_bessel_i check_cyl_bessel_i.cc

check_cyl_bessel_j: check_cyl_bessel_j.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_bessel_j check_cyl_bessel_j.cc

check_cyl_bessel_j_asymp: check_cyl_bessel_j_asymp.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_bessel_j_asymp check_cyl_bessel_j_asymp.cc

check_cyl_bessel_k: check_cyl_bessel_k.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_bessel_k check_cyl_bessel_k.cc

check_cyl_neumann: check_cyl_neumann.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_neumann check_cyl_neumann.cc

check_cyl_neumann_asymp: check_cyl_neumann_asymp.cc
	$(CXX) -D__TEST_DEBUG -o check_cyl_neumann_asymp check_cyl_neumann_asymp.cc

check_ellint_1: check_ellint_1.cc
	$(CXX) -D__TEST_DEBUG -o check_ellint_1 check_ellint_1.cc

check_ellint_2: check_ellint_2.cc
	$(CXX) -D__TEST_DEBUG -o check_ellint_2 check_ellint_2.cc

check_ellint_3: check_ellint_3.cc
	$(CXX) -D__TEST_DEBUG -o check_ellint_3 check_ellint_3.cc

check_expint_neg: check_expint_neg.cc
	$(CXX) -D__TEST_DEBUG -o check_expint_neg check_expint_neg.cc

check_expint_pos: check_expint_pos.cc
	$(CXX) -D__TEST_DEBUG -o check_expint_pos check_expint_pos.cc

check_hyperg: check_hyperg.cc
	$(CXX) -D__TEST_DEBUG -o check_hyperg check_hyperg.cc

check_laguerre: check_laguerre.cc
	$(CXX) -D__TEST_DEBUG -o check_laguerre check_laguerre.cc

check_legendre: check_legendre.cc
	$(CXX) -D__TEST_DEBUG -o check_legendre check_legendre.cc

check_riemann_zeta_neg: check_riemann_zeta_neg.cc
	$(CXX) -D__TEST_DEBUG -o check_riemann_zeta_neg check_riemann_zeta_neg.cc

check_riemann_zeta_pos: check_riemann_zeta_pos.cc
	$(CXX) -D__TEST_DEBUG -o check_riemann_zeta_pos check_riemann_zeta_pos.cc

check_sph_bessel: check_sph_bessel.cc
	$(CXX) -D__TEST_DEBUG -o check_sph_bessel check_sph_bessel.cc

check_sph_legendre: check_sph_legendre.cc
	$(CXX) -D__TEST_DEBUG -o check_sph_legendre check_sph_legendre.cc

check_sph_neumann: check_sph_neumann.cc
	$(CXX) -D__TEST_DEBUG -o check_sph_neumann check_sph_neumann.cc


tarball:
	mkdir tr29124_test
	cp Makefile cmath *.* tr29124_test
	tar -cvf tr29124_test.tar tr29124_test
	bzip2 -f tr29124_test.tar
	md5sum tr29124_test.tar.bz2 > tr29124_test.tar.bz2.md5
	rm -rf tr29124_test

clean:
	rm -f tr29124_*_[fdl].txt
	rm -f gsl_*_[fdl].txt
	rm -f diff_*_[fdl].txt
	rm -f $(BINS)
	rm -f tr29124_test.tar.bz2*

