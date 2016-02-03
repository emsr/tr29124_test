
# -Wconversion

SUFFIX = _specfun

CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
CXX_SRC_DIR = $(HOME)/gcc$(SUFFIX)

CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/6.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = $(CXX_SRC_DIR)/libstdc++-v3/testsuite/util

GSL_DIR = /usr/local
GSL_INC_DIR = $(GSL_DIR)/include
GSL_LIB_DIR = $(GSL_DIR)/lib

GSL_EXT_DIR = $(HOME)/tr29124_test/gslextras
GSL_FRESNEL_DIR = $(GSL_EXT_DIR)/Fresnel

GSL_LIBS = $(GSL_FRESNEL_DIR)/fresnel.c -L$(GSL_LIB_DIR) -lgsl -lgslcblas -ljacobi

BOOST_DIR = /usr/local
BOOST_INC_DIR = $(BOOST_DIR)/include
BOOST_LIB_DIR = $(BOOST_DIR)/lib

#BOOST_LIBS = -L$(BOOST_LIB_DIR) -lboost_math_tools -lboost_math_tr1f -lboost_math_tr1l -lboost_math_tr1

BINS = test_special_function \
       diff_special_function \
       testcase_old \
       testcase \
       airy_toy \
       hankel_toy \
       hankel_toy128 \
       hankel_toy_new \
       test_airy \
       test_csint \
       test_fresnel \
       test_hermite \
       test_limits \
       test_cmath \
       test_bessel \
       test_nric_bessel \
       test_legendre \
       diff_local_special_function \
       test_local_special_function

all: diff_special_function \
     test_special_function \
     testcase_old \
     testcase \
     airy_toy \
     hankel_toy \
     hankel_toy128 \
     hankel_toy_new \
     test_airy \
     test_csint \
     test_fresnel \
     test_hermite \
     test_limits \
     test_cmath \
     test_bessel \
     test_nric_bessel \
     test_legendre \
     diff_local_special_function \
     test_local_special_function

testcases_old: testcase_old
	LD_LIBRARY_PATH=/home/ed/bin_specfun/lib64:$(GSL_LIB_DIR):$$LD_LIBRARY_PATH ./testcase_old

testcases: testcase
	LD_LIBRARY_PATH=/home/ed/bin_specfun/lib64:$(GSL_LIB_DIR):$$LD_LIBRARY_PATH ./testcase

diffs: diff_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./diff_special_function > diff_special_function.txt

tests: test_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_special_function > test_special_function.txt

test:
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_airy > test_airy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_csint > test_csint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_fresnel > test_fresnel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hermite > test_hermite.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_limits > test_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_cmath > test_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bessel > test_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_nric_bessel > test_nric_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./airy_toy > airy_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy > hankel_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy128 > hankel_toy128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy_new > hankel_toy_new.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_legendre > test_legendre.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./diff_local_special_function > diff_local_special_function.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_local_special_function > test_local_special_function.txt

check: \
       check_airy_ai \
       check_airy_bi \
       check_assoc_laguerre \
       check_assoc_legendre \
       check_beta \
       check_bincoef \
       check_comp_ellint_1 \
       check_comp_ellint_2 \
       check_comp_ellint_3 \
       check_conf_hyperg \
       check_conf_hyperg_lim \
       check_coshint \
       check_cosint \
       check_cyl_bessel_i \
       check_cyl_bessel_j \
       check_cyl_bessel_k \
       check_cyl_hankel_1 \
       check_cyl_hankel_2 \
       check_cyl_neumann \
       check_dawson \
       check_dilog \
       check_ellint_1 \
       check_ellint_2 \
       check_ellint_3 \
       check_ellint_d \
       check_ellint_rc \
       check_ellint_rd \
       check_ellint_rf \
       check_ellint_rj \
       check_expint \
       check_expint_e1 \
       check_factorial \
       check_fresnel_c \
       check_fresnel_s \
       check_gamma_l \
       check_gamma_u \
       check_gegenbauer \
       check_hermite \
       check_heuman_lambda \
       check_hurwitz_zeta \
       check_hyperg \
       check_ibeta \
       check_jacobi \
       check_jacobi_cn \
       check_jacobi_dn \
       check_jacobi_sn \
       check_jacobi_zeta \
       check_laguerre \
       check_lbincoef \
       check_ldouble_factorial \
       check_legendre \
       check_legendre \
       check_legendre_q \
       check_lfactorial \
       check_lpochhammer_l \
       check_lpochhammer_u \
       check_pochhammer_l \
       check_pochhammer_u \
       check_psi \
       check_radpoly \
       check_riemann_zeta \
       check_sinc \
       check_sinc_pi \
       check_sinhint \
       check_sinint \
       check_sph_bessel \
       check_sph_bessel_i \
       check_sph_bessel_k \
       check_sph_hankel_1 \
       check_sph_hankel_2 \
       check_sph_harmonic \
       check_sph_legendre \
       check_sph_neumann \
       check_zernike


test_special_function: test_special_function.cpp gsl_wrap.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o test_special_function -I$(GSL_INC_DIR) test_special_function.cpp gsl_wrap.cpp -lquadmath $(GSL_LIBS)

diff_special_function: diff_special_function.cpp gsl_wrap.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o diff_special_function -I$(GSL_INC_DIR) diff_special_function.cpp gsl_wrap.cpp -lquadmath $(GSL_LIBS)

#  You need gnu to get __float128!
test_local_special_function: test_special_function.cpp gsl_wrap.cpp test_func.tcc sf_*.tcc
	$(HOME)/bin/bin/g++ -std=gnu++14 -g -DLOCAL -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -I$(HOME)/gcc_specfun/libstdc++-v3/include -I$(GSL_INC_DIR) -o test_local_special_function test_special_function.cpp gsl_wrap.cpp $(GSL_LIBS) -lquadmath

diff_local_special_function: diff_special_function.cpp gsl_wrap.cpp test_func.tcc sf_*.tcc
	$(HOME)/bin/bin/g++ -std=gnu++14 -g -DLOCAL -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -I$(HOME)/gcc_specfun/libstdc++-v3/include -I$(GSL_INC_DIR) -o diff_local_special_function diff_special_function.cpp gsl_wrap.cpp $(GSL_LIBS) -lquadmath

testcase_old: testcase_old.cpp testcase_old.tcc gsl_wrap.h gsl_wrap.cpp $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o testcase_old -I$(GSL_INC_DIR) testcase_old.cpp gsl_wrap.cpp $(GSL_LIBS)

testcase: testcase.cpp testcase.tcc gsl_wrap.h gsl_wrap.cpp boost_wrap.h boost_wrap.cpp $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o testcase -I$(GSL_INC_DIR) -I$(BOOST_INC_DIR) testcase.cpp gsl_wrap.cpp boost_wrap.cpp $(GSL_LIBS) $(BOOST_LIBS)

test_limits: test_limits.cpp
	$(CXX) -o test_limits test_limits.cpp

test_cmath: test_cmath.cpp
	$(CXX) -o test_cmath test_cmath.cpp

test_airy: test_airy.cpp sf_airy.tcc gsl_wrap.cpp
	$(CXX) -o test_airy -I$(GSL_INC_DIR) test_airy.cpp gsl_wrap.cpp $(GSL_LIBS)

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


airy_toy: airy_toy.cpp
	$(CXX) -o airy_toy airy_toy.cpp -lquadmath

hankel_toy: hankel_toy.cpp
	$(CXX) -o hankel_toy hankel_toy.cpp

hankel_toy128: hankel_toy128.cpp
	$(CXX) -o hankel_toy128 hankel_toy128.cpp -lquadmath

hankel_toy_new: hankel_toy_new.cpp
	$(CXX) -o hankel_toy_new hankel_toy_new.cpp -lquadmath


check_airy_ai: check_airy_ai.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_airy_ai check_airy_ai.cc

check_airy_bi: check_airy_bi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_airy_bi check_airy_bi.cc

check_assoc_laguerre: check_assoc_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_assoc_laguerre check_assoc_laguerre.cc

check_assoc_legendre: check_assoc_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_assoc_legendre check_assoc_legendre.cc

check_beta: check_beta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_beta check_beta.cc

check_bincoef: check_bincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_bincoef check_bincoef.cc

check_comp_ellint_1: check_comp_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_comp_ellint_1 check_comp_ellint_1.cc

check_comp_ellint_2: check_comp_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_comp_ellint_2 check_comp_ellint_2.cc

check_comp_ellint_3: check_comp_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_comp_ellint_3 check_comp_ellint_3.cc

check_conf_hyperg: check_conf_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_conf_hyperg check_conf_hyperg.cc

check_conf_hyperg_lim: check_conf_hyperg_lim.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_conf_hyperg_lim check_conf_hyperg_lim.cc

check_coshint: check_coshint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_coshint check_coshint.cc

check_cosint: check_cosint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cosint check_cosint.cc

check_cyl_bessel_i: check_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_bessel_i check_cyl_bessel_i.cc

check_cyl_bessel_j: check_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_bessel_j check_cyl_bessel_j.cc

check_cyl_bessel_k: check_cyl_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_bessel_k check_cyl_bessel_k.cc

check_cyl_hankel_1: check_cyl_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_hankel_1 check_cyl_hankel_1.cc

check_cyl_hankel_2: check_cyl_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_hankel_2 check_cyl_hankel_2.cc

check_cyl_neumann: check_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_cyl_neumann check_cyl_neumann.cc

check_dawson: check_dawson.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_dawson check_dawson.cc

check_dilog: check_dilog.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_dilog check_dilog.cc

check_ellint_1: check_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_1 check_ellint_1.cc

check_ellint_2: check_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_2 check_ellint_2.cc

check_ellint_3: check_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_3 check_ellint_3.cc

check_ellint_d: check_ellint_d.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_d check_ellint_d.cc

check_ellint_rc: check_ellint_rc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_rc check_ellint_rc.cc

check_ellint_rd: check_ellint_rd.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_rd check_ellint_rd.cc

check_ellint_rf: check_ellint_rf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_rf check_ellint_rf.cc

check_ellint_rj: check_ellint_rj.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ellint_rj check_ellint_rj.cc

check_expint: check_expint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_expint check_expint.cc

check_expint_e1: check_expint_e1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_expint_e1 check_expint_e1.cc

check_factorial: check_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_factorial check_factorial.cc

check_fresnel_c: check_fresnel_c.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_fresnel_c check_fresnel_c.cc

check_fresnel_s: check_fresnel_s.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_fresnel_s check_fresnel_s.cc

check_gamma_l: check_gamma_l.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_gamma_l check_gamma_l.cc

check_gamma_u: check_gamma_u.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_gamma_u check_gamma_u.cc

check_gegenbauer: check_gegenbauer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_gegenbauer check_gegenbauer.cc

check_hermite: check_hermite.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_hermite check_hermite.cc

check_heuman_lambda: check_heuman_lambda.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_heuman_lambda check_heuman_lambda.cc

check_hurwitz_zeta: check_hurwitz_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_hurwitz_zeta check_hurwitz_zeta.cc

check_hyperg: check_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_hyperg check_hyperg.cc

check_ibeta: check_ibeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ibeta check_ibeta.cc

check_jacobi: check_jacobi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_jacobi check_jacobi.cc

check_jacobi_cn: check_jacobi_cn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_jacobi_cn check_jacobi_cn.cc

check_jacobi_dn: check_jacobi_dn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_jacobi_dn check_jacobi_dn.cc

check_jacobi_sn: check_jacobi_sn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_jacobi_sn check_jacobi_sn.cc

check_jacobi_zeta: check_jacobi_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_jacobi_zeta check_jacobi_zeta.cc

check_laguerre: check_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_laguerre check_laguerre.cc

check_lbincoef: check_lbincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_lbincoef check_lbincoef.cc

check_ldouble_factorial: check_ldouble_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_ldouble_factorial check_ldouble_factorial.cc

check_legendre: check_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_legendre check_legendre.cc

check_legendre_q: check_legendre_q.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_legendre_q check_legendre_q.cc

check_lfactorial: check_lfactorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_lfactorial check_lfactorial.cc

check_lpochhammer_l: check_lpochhammer_l.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_lpochhammer_l check_lpochhammer_l.cc

check_lpochhammer_u: check_lpochhammer_u.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_lpochhammer_u check_lpochhammer_u.cc

check_pochhammer_l: check_pochhammer_l.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_pochhammer_l check_pochhammer_l.cc

check_pochhammer_u: check_pochhammer_u.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_pochhammer_u check_pochhammer_u.cc

check_psi: check_psi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_psi check_psi.cc

check_radpoly: check_radpoly.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_radpoly check_radpoly.cc

check_riemann_zeta: check_riemann_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_riemann_zeta check_riemann_zeta.cc

check_sinc: check_sinc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sinc check_sinc.cc

check_sinc_pi: check_sinc_pi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sinc_pi check_sinc_pi.cc

check_sinhint: check_sinhint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sinhint check_sinhint.cc

check_sinint: check_sinint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sinint check_sinint.cc

check_sph_bessel: check_sph_bessel.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_bessel check_sph_bessel.cc

check_sph_bessel_i: check_sph_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_bessel_i check_sph_bessel_i.cc

check_sph_bessel_k: check_sph_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_bessel_k check_sph_bessel_k.cc

check_sph_hankel_1: check_sph_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_hankel_1 check_sph_hankel_1.cc

check_sph_hankel_2: check_sph_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_hankel_2 check_sph_hankel_2.cc

check_sph_harmonic: check_sph_harmonic.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_harmonic check_sph_harmonic.cc

check_sph_legendre: check_sph_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_legendre check_sph_legendre.cc

check_sph_neumann: check_sph_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_sph_neumann check_sph_neumann.cc

check_zernike: check_zernike.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D__TEST_DEBUG -o check_zernike check_zernike.cc


tarball:
	mkdir tr29124
	cp Makefile cmath *.* tr29124
	tar -cvf tr29124.tar tr29124
	bzip2 -f tr29124.tar
	md5sum tr29124.tar.bz2 > tr29124.tar.bz2.md5
	rm -rf tr29124

clean:
	rm -f std_*_[fdl].txt
	rm -f gsl_*_[fdl].txt
	rm -f diff_*_[fdl].txt
	rm -f $(BINS)
	rm -f tr29124.tar.bz2*

