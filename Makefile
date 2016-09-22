
# -Wconversion

SUFFIX = _tr29124

CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
CXX_SRC_DIR = $(HOME)/gcc$(SUFFIX)

CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-compare-reals
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-compare-reals
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/7.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = $(CXX_SRC_DIR)/libstdc++-v3/testsuite/util

INC_DIR = bits

LERCH_DIR = lerchphi/Source

FAD_DIR = Faddeeva

GSL_DIR = /usr/local
GSL_INC_DIR = $(GSL_DIR)/include
GSL_LIB_DIR = $(GSL_DIR)/lib

GSL_EXT_DIR = $(HOME)/tr29124_test/gslextras
GSL_FRESNEL_DIR = $(GSL_EXT_DIR)/Fresnel
GSL_JACOBI_DIR = $(GSL_EXT_DIR)/Jacobi/jacobi-0.9.2/src
GSL_HERMITE_DIR = $(GSL_EXT_DIR)/Hermite

GSL_LIBS = $(GSL_FRESNEL_DIR)/fresnel.c $(GSL_JACOBI_DIR)/jacobi.c $(GSL_HERMITE_DIR)/gsl_sf_hermite.c -L$(GSL_LIB_DIR) -lgsl -lgslcblas

BOOST_DIR = /usr/local
BOOST_INC_DIR = $(BOOST_DIR)/include
BOOST_LIB_DIR = $(BOOST_DIR)/lib

CHECK_DIR = $(HOME)/tr29124_test/check

#BOOST_LIBS = -L$(BOOST_LIB_DIR) -lboost_math_tools -lboost_math_tr1f -lboost_math_tr1l -lboost_math_tr1

BINS = diff_special_function \
       test_special_function \
       testcase2 \
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

CHECKS = ${CHECK_DIR}/check_airy_ai \
	 ${CHECK_DIR}/check_airy_bi \
	 ${CHECK_DIR}/check_assoc_laguerre \
	 ${CHECK_DIR}/check_assoc_legendre \
	 ${CHECK_DIR}/check_beta \
	 ${CHECK_DIR}/check_bincoef \
	 ${CHECK_DIR}/check_chi \
	 ${CHECK_DIR}/check_comp_ellint_1 \
	 ${CHECK_DIR}/check_comp_ellint_2 \
	 ${CHECK_DIR}/check_comp_ellint_3 \
	 ${CHECK_DIR}/check_comp_ellint_d \
	 ${CHECK_DIR}/check_conf_hyperg \
	 ${CHECK_DIR}/check_conf_hyperg_lim \
	 ${CHECK_DIR}/check_coshint \
	 ${CHECK_DIR}/check_cosint \
	 ${CHECK_DIR}/check_cyl_bessel_i \
	 ${CHECK_DIR}/check_cyl_bessel_j \
	 ${CHECK_DIR}/check_cyl_bessel_k \
	 ${CHECK_DIR}/check_cyl_hankel_1 \
	 ${CHECK_DIR}/check_cyl_hankel_2 \
	 ${CHECK_DIR}/check_cyl_neumann \
	 ${CHECK_DIR}/check_dawson \
	 ${CHECK_DIR}/check_dilog \
	 ${CHECK_DIR}/check_ellint_1 \
	 ${CHECK_DIR}/check_ellint_2 \
	 ${CHECK_DIR}/check_ellint_3 \
	 ${CHECK_DIR}/check_ellint_d \
	 ${CHECK_DIR}/check_ellint_rc \
	 ${CHECK_DIR}/check_ellint_rd \
	 ${CHECK_DIR}/check_ellint_rf \
	 ${CHECK_DIR}/check_ellint_rg \
	 ${CHECK_DIR}/check_ellint_rj \
	 ${CHECK_DIR}/check_expint \
	 ${CHECK_DIR}/check_expint_en \
	 ${CHECK_DIR}/check_factorial \
	 ${CHECK_DIR}/check_fresnel_c \
	 ${CHECK_DIR}/check_fresnel_s \
	 ${CHECK_DIR}/check_gegenbauer \
	 ${CHECK_DIR}/check_hermite \
	 ${CHECK_DIR}/check_heuman_lambda \
	 ${CHECK_DIR}/check_hurwitz_zeta \
	 ${CHECK_DIR}/check_hyperg \
	 ${CHECK_DIR}/check_ibeta \
	 ${CHECK_DIR}/check_jacobi \
	 ${CHECK_DIR}/check_jacobi_cn \
	 ${CHECK_DIR}/check_jacobi_dn \
	 ${CHECK_DIR}/check_jacobi_sn \
	 ${CHECK_DIR}/check_jacobi_zeta \
	 ${CHECK_DIR}/check_laguerre \
	 ${CHECK_DIR}/check_lbincoef \
	 ${CHECK_DIR}/check_ldouble_factorial \
	 ${CHECK_DIR}/check_legendre \
	 ${CHECK_DIR}/check_legendre \
	 ${CHECK_DIR}/check_legendre_q \
	 ${CHECK_DIR}/check_lfactorial \
	 ${CHECK_DIR}/check_lpochhammer_lower \
	 ${CHECK_DIR}/check_lpochhammer \
	 ${CHECK_DIR}/check_owens_t \
	 ${CHECK_DIR}/check_pgamma \
	 ${CHECK_DIR}/check_pochhammer_lower \
	 ${CHECK_DIR}/check_pochhammer \
	 ${CHECK_DIR}/check_psi \
	 ${CHECK_DIR}/check_qgamma \
	 ${CHECK_DIR}/check_radpoly \
	 ${CHECK_DIR}/check_riemann_zeta \
	 ${CHECK_DIR}/check_shi \
	 ${CHECK_DIR}/check_sinc \
	 ${CHECK_DIR}/check_sinc_pi \
	 ${CHECK_DIR}/check_sinhint \
	 ${CHECK_DIR}/check_sinint \
	 ${CHECK_DIR}/check_sph_bessel \
	 ${CHECK_DIR}/check_sph_bessel_i \
	 ${CHECK_DIR}/check_sph_bessel_k \
	 ${CHECK_DIR}/check_sph_hankel_1 \
	 ${CHECK_DIR}/check_sph_hankel_2 \
	 ${CHECK_DIR}/check_sph_harmonic \
	 ${CHECK_DIR}/check_sph_legendre \
	 ${CHECK_DIR}/check_sph_neumann \
	 ${CHECK_DIR}/check_tgamma_lower \
	 ${CHECK_DIR}/check_tgamma \
	 ${CHECK_DIR}/check_theta_1 \
	 ${CHECK_DIR}/check_theta_2 \
	 ${CHECK_DIR}/check_theta_3 \
	 ${CHECK_DIR}/check_theta_4 \
	 ${CHECK_DIR}/check_zernike \
	 ${CHECK_DIR}/complex_ellint_rc \
	 ${CHECK_DIR}/complex_ellint_rd \
	 ${CHECK_DIR}/complex_ellint_rf \
	 ${CHECK_DIR}/complex_ellint_rg \
	 ${CHECK_DIR}/complex_ellint_rj \
	 ${CHECK_DIR}/check_clausen_c \
	 ${CHECK_DIR}/pr56216_cyl_hankel_1 \
	 ${CHECK_DIR}/pr56216_cyl_hankel_2 \
	 ${CHECK_DIR}/pr56216_cyl_bessel_i \
	 ${CHECK_DIR}/origin_bessel_j \
	 ${CHECK_DIR}/origin_cyl_neumann

all: diff_special_function \
     test_special_function \
     testcase2 \
     testcase \
     airy_toy \
     hankel_toy \
     hankel_toy128 \
     hankel_toy_new \
     test_airy \
     test_csint \
     test_Faddeeva \
     test_fresnel \
     test_hermite \
     test_limits \
     test_cmath \
     test_bessel \
     test_nric_bessel \
     test_legendre \
     diff_local_special_function \
     test_local_special_function


docs: bits/*
	rm -rf html/*
	rm -rf latex/*
	doxygen
	cd latex && make

testcases2: testcase2
	LD_LIBRARY_PATH=/home/ed/bin$(SUFFIX)/lib64:$(GSL_LIB_DIR):$$LD_LIBRARY_PATH ./testcase2

testcases: testcase
	LD_LIBRARY_PATH=/home/ed/bin$(SUFFIX)/lib64:$(GSL_LIB_DIR):$$LD_LIBRARY_PATH ./testcase

diffs: diff_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./diff_special_function > diff_special_function.txt

tests: test_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_special_function > test_special_function.txt

test:
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_airy > test_airy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_csint > test_csint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(FAD_DIR)/test_Faddeeva > $(FAD_DIR)/test_Faddeeva.txt
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

check: $(CHECKS)
	echo "Beginning executions of checks..." > $(CHECK_DIR)/check_out.txt 2> $(CHECK_DIR)/check_err.txt
	echo "check_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bincoef" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_bincoef >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_clausen_c" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_clausen_c >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_comp_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_comp_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_comp_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_d" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_comp_ellint_d >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_conf_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_conf_hyperg_lim" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_conf_hyperg_lim >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_coshint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_coshint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cosint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cosint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_k" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dawson" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dawson >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dilog" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dilog >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_d" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_d >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rd" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rd >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rj" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rj >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_expint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint_en" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_expint_en >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_c" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_c >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_s" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_s >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gegenbauer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gegenbauer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hermite" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hermite >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_heuman_lambda" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_heuman_lambda >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hurwitz_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hurwitz_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ibeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ibeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_jacobi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_cn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_jacobi_cn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_dn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_jacobi_dn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_sn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_jacobi_sn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_jacobi_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lbincoef" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lbincoef >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ldouble_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ldouble_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lfactorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lfactorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lpochhammer_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lpochhammer_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lpochhammer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lpochhammer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_owens_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_owens_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pochhammer_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pochhammer_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pochhammer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pochhammer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_psi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_psi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_qgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_qgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_radpoly" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_radpoly >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_riemann_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_shi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_shi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinc_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinhint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinhint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_bessel >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel_k" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_harmonic" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_harmonic >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sph_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tgamma_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_tgamma_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_tgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_theta_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_theta_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_theta_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_4" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_theta_4 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_zernike" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_zernike >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_ellint_rc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rd" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_ellint_rd >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_ellint_rf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_ellint_rg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rj" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_ellint_rj >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


test_special_function: test_special_function.cpp wrap_gsl.h wrap_gsl.cpp wrap_boost.h wrap_boost.cpp $(LERCH_DIR)/lerchphi.h $(LERCH_DIR)/lerchphi.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o test_special_function -I$(GSL_INC_DIR) test_special_function.cpp wrap_gsl.cpp wrap_boost.cpp $(LERCH_DIR)/lerchphi.cpp -lquadmath $(GSL_LIBS)

diff_special_function: diff_special_function.cpp wrap_gsl.h wrap_gsl.cpp wrap_boost.h wrap_boost.cpp $(LERCH_DIR)/lerchphi.h $(LERCH_DIR)/lerchphi.cpp test_func.tcc $(CXX_INC_DIR)/sf_*.tcc
	$(CXX) -o diff_special_function -I$(GSL_INC_DIR) diff_special_function.cpp wrap_gsl.cpp wrap_boost.cpp $(LERCH_DIR)/lerchphi.cpp -lquadmath $(GSL_LIBS)

#  You need gnu to get __float128!
test_local_special_function: test_special_function.cpp wrap_gsl.cpp test_func.tcc sf_*.tcc
	$(HOME)/bin/bin/g++ -std=gnu++14 -g -DLOCAL -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -I$(HOME)/gcc$(SUFFIX)/libstdc++-v3/include -I$(GSL_INC_DIR) -o test_local_special_function test_special_function.cpp wrap_gsl.cpp $(GSL_LIBS) -lquadmath

diff_local_special_function: diff_special_function.cpp wrap_gsl.cpp test_func.tcc sf_*.tcc
	$(HOME)/bin/bin/g++ -std=gnu++14 -g -DLOCAL -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -I$(HOME)/gcc$(SUFFIX)/libstdc++-v3/include -I$(GSL_INC_DIR) -o diff_local_special_function diff_special_function.cpp wrap_gsl.cpp $(GSL_LIBS) -lquadmath

testcase2: testcase2.cpp testcase2.tcc wrap_gsl.h wrap_gsl.cpp wrap_boost.h wrap_boost.cpp $(LERCH_DIR)/lerchphi.h $(LERCH_DIR)/lerchphi.cpp wrap_burkhardt.h wrap_burkhardt.cpp burkhardt/special_functions.f90 $(INC_DIR)/sf_*.tcc
	$(CXX17) -o testcase2 -I. -I$(GSL_INC_DIR) -I$(BOOST_INC_DIR) testcase2.cpp wrap_gsl.cpp wrap_boost.cpp wrap_burkhardt.cpp burkhardt/special_functions.f90 $(LERCH_DIR)/lerchphi.cpp $(GSL_LIBS) $(BOOST_LIBS) -lgfortran

testcase: testcase.cpp testcase.tcc wrap_gsl.h wrap_gsl.cpp wrap_boost.h wrap_boost.cpp $(LERCH_DIR)/lerchphi.h $(LERCH_DIR)/lerchphi.cpp wrap_burkhardt.h wrap_burkhardt.cpp burkhardt/special_functions.f90 $(INC_DIR)/sf_*.tcc
	$(CXX17) -o testcase -I. -I$(GSL_INC_DIR) -I$(BOOST_INC_DIR) testcase.cpp wrap_gsl.cpp wrap_boost.cpp wrap_burkhardt.cpp burkhardt/special_functions.f90 $(LERCH_DIR)/lerchphi.cpp $(GSL_LIBS) $(BOOST_LIBS) -lgfortran

test_limits: test_limits.cpp
	$(CXX17) -o test_limits test_limits.cpp

test_cmath: test_cmath.cpp
	$(CXX) -o test_cmath test_cmath.cpp

test_airy: test_airy.cpp sf_airy.tcc wrap_gsl.cpp
	$(CXX) -o test_airy -I$(GSL_INC_DIR) test_airy.cpp wrap_gsl.cpp $(GSL_LIBS)

test_csint: test_csint.cpp csint.tcc
	$(CXX) -o test_csint test_csint.cpp

test_Faddeeva: $(FAD_DIR)/Faddeeva.h $(FAD_DIR)/Faddeeva.cpp
	$(CXX) -DTEST_FADDEEVA -o $(FAD_DIR)/test_Faddeeva  $(FAD_DIR)/Faddeeva.cpp

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


${CHECK_DIR}/check_airy_ai: ${CHECK_DIR}/check_airy_ai.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_ai ${CHECK_DIR}/check_airy_ai.cc

${CHECK_DIR}/check_airy_bi: ${CHECK_DIR}/check_airy_bi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_bi ${CHECK_DIR}/check_airy_bi.cc

${CHECK_DIR}/check_assoc_laguerre: ${CHECK_DIR}/check_assoc_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_laguerre ${CHECK_DIR}/check_assoc_laguerre.cc

${CHECK_DIR}/check_assoc_legendre: ${CHECK_DIR}/check_assoc_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_legendre ${CHECK_DIR}/check_assoc_legendre.cc

${CHECK_DIR}/check_beta: ${CHECK_DIR}/check_beta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_beta ${CHECK_DIR}/check_beta.cc

${CHECK_DIR}/check_bincoef: ${CHECK_DIR}/check_bincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_bincoef ${CHECK_DIR}/check_bincoef.cc

${CHECK_DIR}/check_chi: ${CHECK_DIR}/check_chi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chi ${CHECK_DIR}/check_chi.cc

${CHECK_DIR}/check_clausen_c: ${CHECK_DIR}/check_clausen_c.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_clausen_c ${CHECK_DIR}/check_clausen_c.cc

${CHECK_DIR}/check_comp_ellint_1: ${CHECK_DIR}/check_comp_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_1 ${CHECK_DIR}/check_comp_ellint_1.cc

${CHECK_DIR}/check_comp_ellint_2: ${CHECK_DIR}/check_comp_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_2 ${CHECK_DIR}/check_comp_ellint_2.cc

${CHECK_DIR}/check_comp_ellint_3: ${CHECK_DIR}/check_comp_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_3 ${CHECK_DIR}/check_comp_ellint_3.cc

${CHECK_DIR}/check_comp_ellint_d: ${CHECK_DIR}/check_comp_ellint_d.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_d ${CHECK_DIR}/check_comp_ellint_d.cc

${CHECK_DIR}/check_conf_hyperg: ${CHECK_DIR}/check_conf_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg ${CHECK_DIR}/check_conf_hyperg.cc

${CHECK_DIR}/check_conf_hyperg_lim: ${CHECK_DIR}/check_conf_hyperg_lim.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg_lim ${CHECK_DIR}/check_conf_hyperg_lim.cc

${CHECK_DIR}/check_coshint: ${CHECK_DIR}/check_coshint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_coshint ${CHECK_DIR}/check_coshint.cc

${CHECK_DIR}/check_cosint: ${CHECK_DIR}/check_cosint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cosint ${CHECK_DIR}/check_cosint.cc

${CHECK_DIR}/check_cyl_bessel_i: ${CHECK_DIR}/check_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_i ${CHECK_DIR}/check_cyl_bessel_i.cc

${CHECK_DIR}/check_cyl_bessel_j: ${CHECK_DIR}/check_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_j ${CHECK_DIR}/check_cyl_bessel_j.cc

${CHECK_DIR}/check_cyl_bessel_k: ${CHECK_DIR}/check_cyl_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_k ${CHECK_DIR}/check_cyl_bessel_k.cc

${CHECK_DIR}/check_cyl_hankel_1: ${CHECK_DIR}/check_cyl_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_1 ${CHECK_DIR}/check_cyl_hankel_1.cc

${CHECK_DIR}/check_cyl_hankel_2: ${CHECK_DIR}/check_cyl_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_2 ${CHECK_DIR}/check_cyl_hankel_2.cc

${CHECK_DIR}/check_cyl_neumann: ${CHECK_DIR}/check_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_neumann ${CHECK_DIR}/check_cyl_neumann.cc

${CHECK_DIR}/check_dawson: ${CHECK_DIR}/check_dawson.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dawson ${CHECK_DIR}/check_dawson.cc

${CHECK_DIR}/check_dilog: ${CHECK_DIR}/check_dilog.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dilog ${CHECK_DIR}/check_dilog.cc

${CHECK_DIR}/check_ellint_1: ${CHECK_DIR}/check_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_1 ${CHECK_DIR}/check_ellint_1.cc

${CHECK_DIR}/check_ellint_2: ${CHECK_DIR}/check_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_2 ${CHECK_DIR}/check_ellint_2.cc

${CHECK_DIR}/check_ellint_3: ${CHECK_DIR}/check_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_3 ${CHECK_DIR}/check_ellint_3.cc

${CHECK_DIR}/check_ellint_d: ${CHECK_DIR}/check_ellint_d.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_d ${CHECK_DIR}/check_ellint_d.cc

${CHECK_DIR}/check_ellint_rc: ${CHECK_DIR}/check_ellint_rc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rc ${CHECK_DIR}/check_ellint_rc.cc

${CHECK_DIR}/check_ellint_rd: ${CHECK_DIR}/check_ellint_rd.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rd ${CHECK_DIR}/check_ellint_rd.cc

${CHECK_DIR}/check_ellint_rf: ${CHECK_DIR}/check_ellint_rf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rf ${CHECK_DIR}/check_ellint_rf.cc

${CHECK_DIR}/check_ellint_rg: ${CHECK_DIR}/check_ellint_rg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rg ${CHECK_DIR}/check_ellint_rg.cc

${CHECK_DIR}/check_ellint_rj: ${CHECK_DIR}/check_ellint_rj.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rj ${CHECK_DIR}/check_ellint_rj.cc

${CHECK_DIR}/check_expint: ${CHECK_DIR}/check_expint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint ${CHECK_DIR}/check_expint.cc

${CHECK_DIR}/check_expint_en: ${CHECK_DIR}/check_expint_en.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint_en ${CHECK_DIR}/check_expint_en.cc

${CHECK_DIR}/check_factorial: ${CHECK_DIR}/check_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_factorial ${CHECK_DIR}/check_factorial.cc

${CHECK_DIR}/check_fresnel_c: ${CHECK_DIR}/check_fresnel_c.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_c ${CHECK_DIR}/check_fresnel_c.cc

${CHECK_DIR}/check_fresnel_s: ${CHECK_DIR}/check_fresnel_s.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_s ${CHECK_DIR}/check_fresnel_s.cc

${CHECK_DIR}/check_gegenbauer: ${CHECK_DIR}/check_gegenbauer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gegenbauer ${CHECK_DIR}/check_gegenbauer.cc

${CHECK_DIR}/check_hermite: ${CHECK_DIR}/check_hermite.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hermite ${CHECK_DIR}/check_hermite.cc

${CHECK_DIR}/check_heuman_lambda: ${CHECK_DIR}/check_heuman_lambda.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_heuman_lambda ${CHECK_DIR}/check_heuman_lambda.cc

${CHECK_DIR}/check_hurwitz_zeta: ${CHECK_DIR}/check_hurwitz_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hurwitz_zeta ${CHECK_DIR}/check_hurwitz_zeta.cc

${CHECK_DIR}/check_hyperg: ${CHECK_DIR}/check_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hyperg ${CHECK_DIR}/check_hyperg.cc

${CHECK_DIR}/check_ibeta: ${CHECK_DIR}/check_ibeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ibeta ${CHECK_DIR}/check_ibeta.cc

${CHECK_DIR}/check_jacobi: ${CHECK_DIR}/check_jacobi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi ${CHECK_DIR}/check_jacobi.cc

${CHECK_DIR}/check_jacobi_cn: ${CHECK_DIR}/check_jacobi_cn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_cn ${CHECK_DIR}/check_jacobi_cn.cc

${CHECK_DIR}/check_jacobi_dn: ${CHECK_DIR}/check_jacobi_dn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_dn ${CHECK_DIR}/check_jacobi_dn.cc

${CHECK_DIR}/check_jacobi_sn: ${CHECK_DIR}/check_jacobi_sn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_sn ${CHECK_DIR}/check_jacobi_sn.cc

${CHECK_DIR}/check_jacobi_zeta: ${CHECK_DIR}/check_jacobi_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_zeta ${CHECK_DIR}/check_jacobi_zeta.cc

${CHECK_DIR}/check_laguerre: ${CHECK_DIR}/check_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_laguerre ${CHECK_DIR}/check_laguerre.cc

${CHECK_DIR}/check_lbincoef: ${CHECK_DIR}/check_lbincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lbincoef ${CHECK_DIR}/check_lbincoef.cc

${CHECK_DIR}/check_ldouble_factorial: ${CHECK_DIR}/check_ldouble_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ldouble_factorial ${CHECK_DIR}/check_ldouble_factorial.cc

${CHECK_DIR}/check_legendre: ${CHECK_DIR}/check_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre ${CHECK_DIR}/check_legendre.cc

${CHECK_DIR}/check_legendre_q: ${CHECK_DIR}/check_legendre_q.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre_q ${CHECK_DIR}/check_legendre_q.cc

${CHECK_DIR}/check_lfactorial: ${CHECK_DIR}/check_lfactorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfactorial ${CHECK_DIR}/check_lfactorial.cc

${CHECK_DIR}/check_lpochhammer_lower: ${CHECK_DIR}/check_lpochhammer_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lpochhammer_lower ${CHECK_DIR}/check_lpochhammer_lower.cc

${CHECK_DIR}/check_lpochhammer: ${CHECK_DIR}/check_lpochhammer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lpochhammer ${CHECK_DIR}/check_lpochhammer.cc

${CHECK_DIR}/check_owens_t: ${CHECK_DIR}/check_owens_t.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_owens_t ${CHECK_DIR}/check_owens_t.cc

${CHECK_DIR}/check_pgamma: ${CHECK_DIR}/check_pgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pgamma ${CHECK_DIR}/check_pgamma.cc

${CHECK_DIR}/check_pochhammer_lower: ${CHECK_DIR}/check_pochhammer_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pochhammer_lower ${CHECK_DIR}/check_pochhammer_lower.cc

${CHECK_DIR}/check_pochhammer: ${CHECK_DIR}/check_pochhammer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pochhammer ${CHECK_DIR}/check_pochhammer.cc

${CHECK_DIR}/check_psi: ${CHECK_DIR}/check_psi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_psi ${CHECK_DIR}/check_psi.cc

${CHECK_DIR}/check_qgamma: ${CHECK_DIR}/check_qgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_qgamma ${CHECK_DIR}/check_qgamma.cc

${CHECK_DIR}/check_radpoly: ${CHECK_DIR}/check_radpoly.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_radpoly ${CHECK_DIR}/check_radpoly.cc

${CHECK_DIR}/check_riemann_zeta: ${CHECK_DIR}/check_riemann_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_riemann_zeta ${CHECK_DIR}/check_riemann_zeta.cc

${CHECK_DIR}/check_shi: ${CHECK_DIR}/check_shi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_shi ${CHECK_DIR}/check_shi.cc

${CHECK_DIR}/check_sinc: ${CHECK_DIR}/check_sinc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc ${CHECK_DIR}/check_sinc.cc

${CHECK_DIR}/check_sinc_pi: ${CHECK_DIR}/check_sinc_pi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc_pi ${CHECK_DIR}/check_sinc_pi.cc

${CHECK_DIR}/check_sinhint: ${CHECK_DIR}/check_sinhint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinhint ${CHECK_DIR}/check_sinhint.cc

${CHECK_DIR}/check_sinint: ${CHECK_DIR}/check_sinint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinint ${CHECK_DIR}/check_sinint.cc

${CHECK_DIR}/check_sph_bessel: ${CHECK_DIR}/check_sph_bessel.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel ${CHECK_DIR}/check_sph_bessel.cc

${CHECK_DIR}/check_sph_bessel_i: ${CHECK_DIR}/check_sph_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_i ${CHECK_DIR}/check_sph_bessel_i.cc

${CHECK_DIR}/check_sph_bessel_k: ${CHECK_DIR}/check_sph_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_k ${CHECK_DIR}/check_sph_bessel_k.cc

${CHECK_DIR}/check_sph_hankel_1: ${CHECK_DIR}/check_sph_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_1 ${CHECK_DIR}/check_sph_hankel_1.cc

${CHECK_DIR}/check_sph_hankel_2: ${CHECK_DIR}/check_sph_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_2 ${CHECK_DIR}/check_sph_hankel_2.cc

${CHECK_DIR}/check_sph_harmonic: ${CHECK_DIR}/check_sph_harmonic.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_harmonic ${CHECK_DIR}/check_sph_harmonic.cc

${CHECK_DIR}/check_sph_legendre: ${CHECK_DIR}/check_sph_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_legendre ${CHECK_DIR}/check_sph_legendre.cc

${CHECK_DIR}/check_sph_neumann: ${CHECK_DIR}/check_sph_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_neumann ${CHECK_DIR}/check_sph_neumann.cc

${CHECK_DIR}/check_tgamma_lower: ${CHECK_DIR}/check_tgamma_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma_lower ${CHECK_DIR}/check_tgamma_lower.cc

${CHECK_DIR}/check_tgamma: ${CHECK_DIR}/check_tgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma ${CHECK_DIR}/check_tgamma.cc

${CHECK_DIR}/check_theta_1: ${CHECK_DIR}/check_theta_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_1 ${CHECK_DIR}/check_theta_1.cc

${CHECK_DIR}/check_theta_2: ${CHECK_DIR}/check_theta_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_2 ${CHECK_DIR}/check_theta_2.cc

${CHECK_DIR}/check_theta_3: ${CHECK_DIR}/check_theta_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_3 ${CHECK_DIR}/check_theta_3.cc

${CHECK_DIR}/check_theta_4: ${CHECK_DIR}/check_theta_4.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_4 ${CHECK_DIR}/check_theta_4.cc

${CHECK_DIR}/check_zernike: ${CHECK_DIR}/check_zernike.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_zernike ${CHECK_DIR}/check_zernike.cc

${CHECK_DIR}/complex_ellint_rc: ${CHECK_DIR}/complex_ellint_rc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rc ${CHECK_DIR}/complex_ellint_rc.cc

${CHECK_DIR}/complex_ellint_rd: ${CHECK_DIR}/complex_ellint_rd.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rd ${CHECK_DIR}/complex_ellint_rd.cc

${CHECK_DIR}/complex_ellint_rf: ${CHECK_DIR}/complex_ellint_rf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rf ${CHECK_DIR}/complex_ellint_rf.cc

${CHECK_DIR}/complex_ellint_rg: ${CHECK_DIR}/complex_ellint_rg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rg ${CHECK_DIR}/complex_ellint_rg.cc

${CHECK_DIR}/complex_ellint_rj: ${CHECK_DIR}/complex_ellint_rj.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rj ${CHECK_DIR}/complex_ellint_rj.cc

${CHECK_DIR}/pr56216_cyl_hankel_1: ${CHECK_DIR}/pr56216_cyl_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_1 ${CHECK_DIR}/pr56216_cyl_hankel_1.cc

${CHECK_DIR}/pr56216_cyl_hankel_2: ${CHECK_DIR}/pr56216_cyl_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_2 ${CHECK_DIR}/pr56216_cyl_hankel_2.cc

${CHECK_DIR}/pr56216_cyl_bessel_i: ${CHECK_DIR}/pr56216_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_bessel_i ${CHECK_DIR}/pr56216_cyl_bessel_i.cc

${CHECK_DIR}/origin_bessel_j: ${CHECK_DIR}/origin_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_bessel_j ${CHECK_DIR}/origin_bessel_j.cc

${CHECK_DIR}/origin_cyl_neumann: ${CHECK_DIR}/origin_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_neumann ${CHECK_DIR}/origin_cyl_neumann.cc


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
	rm -f $(CHECKS)
	rm -f tr29124.tar.bz2*

