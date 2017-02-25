
# -Wconversion

SUFFIX = _tr29124

CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
CXX_SRC_DIR = $(HOME)/gcc$(SUFFIX)

GFORTRAN = $(CXX_INST_DIR)/bin/gfortran -g -Wall -Wextra -Wno-compare-reals
GCC = $(CXX_INST_DIR)/bin/gcc -g -Wall -Wextra
CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/7.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
#CXX_TEST_INC_DIR = $(CXX_SRC_DIR)/libstdc++-v3/testsuite/util
CXX_TEST_INC_DIR = .

INC_DIR = bits

CEPHES_DIR = cephes

LERCH_DIR = lerchphi/Source

FAD_DIR = Faddeeva

GSL_DIR = /usr/local
GSL_INC_DIR = $(GSL_DIR)/include
GSL_LIB_DIR = $(GSL_DIR)/lib

GSL_EXT_DIR = $(HOME)/tr29124_test/gslextras
GSL_FRESNEL_DIR = $(GSL_EXT_DIR)/Fresnel
GSL_JACOBI_DIR = $(GSL_EXT_DIR)/Jacobi/jacobi-0.9.2/src
GSL_HERMITE_DIR = $(GSL_EXT_DIR)/Hermite

BOOST_DIR = /usr/local
BOOST_INC_DIR = $(BOOST_DIR)/include
BOOST_LIB_DIR = $(BOOST_DIR)/lib

MATH_DIR = $(HOME)/tr29124_test

CHECK_DIR = $(MATH_DIR)/check

OBJ_DIR = $(MATH_DIR)/obj

LIB_DIR = $(MATH_DIR)

.PHONY: wrappers_debug wrappers_release

wrappers_debug:
	$(MAKE) -C wrappers/debug

wrappers_release:
	$(MAKE) -C wrappers/release

BINS = testcase \
       mpfrcalc \
       diff_special_function \
       test_special_function \
       testcase2 \
       airy_toy \
       hankel_toy \
       hankel_toy128 \
       hankel_toy_new \
       test_airy \
       test_anger_weber \
       test_bernoulli \
       test_bessel \
       test_bessel_asymp \
       test_bessel_iter \
       test_beta \
       test_beta_inc \
       test_binet \
       test_binet_float \
       test_bose_einstein \
       test_chebyshev \
       test_chebyshev_trig \
       test_chebyshev_trig_pi \
       test_clausen \
       test_cmath \
       test_comp_ellint_1 \
       test_complex128 \
       test_complex_gamma \
       test_conf_hyperg_limit \
       test_const \
       test_continued_fraction \
       test_csint \
       test_cyl_hankel \
       test_dawson \
       test_debye \
       test_dilog \
       test_expint \
       test_factorial \
       test_fermi_dirac \
       test_float128 \
       test_fresnel \
       test_gamma \
       test_gamma_ratio \
       test_gamma_reciprocal \
       test_gegenbauer \
       test_hankel \
       test_hankel_real_arg \
       test_hermite \
       test_heuman_lambda \
       test_hurwitz_zeta \
       test_hurwitz_zeta_new \
       test_hyperg \
       test_hypot \
       test_inv_erf \
       test_inv_ibeta \
       test_jacobi \
       test_jacobi_inv \
       test_jacobi_zeta \
       test_kelvin \
       test_laguerre \
       test_legendre \
       test_lentz_continued_fraction \
       test_lerch \
       test_limits \
       test_little_airy \
       test_math_h \
       test_notsospecfun \
       test_nric_bessel \
       test_numeric_limits \
       test_owens_t \
       test_parab_cyl \
       test_pochhammer \
       test_pochhammer_lower \
       test_polylog \
       test_polynomial \
       test_pow_funs \
       test_psi \
       test_rational \
       test_ratpoly \
       test_recursion \
       test_reperiodized_hyper \
       test_reperiodized_trig \
       test_riemann_zeta \
       test_root_finding \
       test_sincos \
       test_sinus_cardinal \
       test_sph_bessel \
       test_sph_hankel \
       test_static_polynomial \
       test_steed_continued_fraction \
       test_struve \
       test_struve_old \
       test_summation \
       test_theta \
       test_trigamma

CHECKS = ${CHECK_DIR}/check_airy_ai \
	 ${CHECK_DIR}/check_airy_bi \
	 ${CHECK_DIR}/check_assoc_laguerre \
	 ${CHECK_DIR}/check_assoc_legendre \
	 ${CHECK_DIR}/check_beta \
	 ${CHECK_DIR}/check_bernoulli \
	 ${CHECK_DIR}/check_bincoef \
	 ${CHECK_DIR}/check_chebyshev_t \
	 ${CHECK_DIR}/check_chebyshev_u \
	 ${CHECK_DIR}/check_chebyshev_v \
	 ${CHECK_DIR}/check_chebyshev_w \
	 ${CHECK_DIR}/check_chi \
	 ${CHECK_DIR}/check_comp_ellint_1 \
	 ${CHECK_DIR}/check_comp_ellint_2 \
	 ${CHECK_DIR}/check_comp_ellint_3 \
	 ${CHECK_DIR}/check_comp_ellint_d \
	 ${CHECK_DIR}/check_conf_hyperg \
	 ${CHECK_DIR}/check_conf_hyperg_lim \
	 ${CHECK_DIR}/check_coshint \
	 ${CHECK_DIR}/check_cosint \
	 ${CHECK_DIR}/check_cos_pi \
	 ${CHECK_DIR}/check_cyl_bessel_i \
	 ${CHECK_DIR}/check_cyl_bessel_j \
	 ${CHECK_DIR}/check_cyl_bessel_k \
	 ${CHECK_DIR}/check_cyl_hankel_1 \
	 ${CHECK_DIR}/check_cyl_hankel_2 \
	 ${CHECK_DIR}/check_cyl_neumann \
	 ${CHECK_DIR}/check_dawson \
	 ${CHECK_DIR}/check_dilog \
	 ${CHECK_DIR}/check_dirichlet_beta \
	 ${CHECK_DIR}/check_dirichlet_eta \
	 ${CHECK_DIR}/check_dirichlet_lambda \
	 ${CHECK_DIR}/check_double_factorial \
	 ${CHECK_DIR}/check_ellint_1 \
	 ${CHECK_DIR}/check_ellint_2 \
	 ${CHECK_DIR}/check_ellint_3 \
	 ${CHECK_DIR}/check_ellint_d \
	 ${CHECK_DIR}/check_ellint_rc \
	 ${CHECK_DIR}/check_ellint_rd \
	 ${CHECK_DIR}/check_ellint_rf \
	 ${CHECK_DIR}/check_ellint_rg \
	 ${CHECK_DIR}/check_ellint_rj \
	 ${CHECK_DIR}/check_ellnome \
	 ${CHECK_DIR}/check_expint \
	 ${CHECK_DIR}/check_expint_en \
	 ${CHECK_DIR}/check_factorial \
	 ${CHECK_DIR}/check_fresnel_c \
	 ${CHECK_DIR}/check_fresnel_s \
	 ${CHECK_DIR}/check_gamma_reciprocal \
	 ${CHECK_DIR}/check_gegenbauer \
	 ${CHECK_DIR}/check_hermite \
	 ${CHECK_DIR}/check_heuman_lambda \
	 ${CHECK_DIR}/check_hurwitz_zeta \
	 ${CHECK_DIR}/check_hyperg \
	 ${CHECK_DIR}/check_ibeta \
	 ${CHECK_DIR}/check_ibetac \
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
	 ${CHECK_DIR}/check_lgamma \
	 ${CHECK_DIR}/check_logistic_cdf \
	 ${CHECK_DIR}/check_logistic_pdf \
	 ${CHECK_DIR}/check_lpochhammer_lower \
	 ${CHECK_DIR}/check_lpochhammer \
	 ${CHECK_DIR}/check_lognormal_cdf \
	 ${CHECK_DIR}/check_lognormal_pdf \
	 ${CHECK_DIR}/check_normal_cdf \
	 ${CHECK_DIR}/check_normal_pdf \
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
	 ${CHECK_DIR}/check_sin_pi \
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
	 ${CHECK_DIR}/complex_airy_ai \
	 ${CHECK_DIR}/complex_airy_bi \
	 ${CHECK_DIR}/deathmatch_comp_ellint \
	 ${CHECK_DIR}/deathmatch_conf_hyperg \
	 ${CHECK_DIR}/deathmatch_conf_hyperg_lim \
	 ${CHECK_DIR}/deathmatch_hyperg \
	 ${CHECK_DIR}/check_clausen_c \
	 ${CHECK_DIR}/pr56216_cyl_hankel_1 \
	 ${CHECK_DIR}/pr56216_cyl_hankel_2 \
	 ${CHECK_DIR}/pr56216_cyl_bessel_i \
	 ${CHECK_DIR}/origin_cyl_bessel_j \
	 ${CHECK_DIR}/origin_cyl_neumann

all: $(BINS)


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

check: $(CHECKS)
	echo "Beginning executions of checks..." > $(CHECK_DIR)/check_out.txt 2> $(CHECK_DIR)/check_err.txt
	echo "check_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bernoulli" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_bernoulli >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bincoef" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_bincoef >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_u" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_u >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_v" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_v >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_w" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_w >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_cos_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cos_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_k" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dawson" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dawson >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dilog" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dilog >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dirichlet_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_eta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dirichlet_eta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_lambda" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_dirichlet_lambda >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_double_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_double_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_d" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_d >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rd" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rd >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rj" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellint_rj >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellnome" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ellnome >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_expint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint_en" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_expint_en >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_c" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_c >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_s" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_s >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_reciprocal" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gamma_reciprocal >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gegenbauer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gegenbauer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hermite" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hermite >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_heuman_lambda" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_heuman_lambda >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hurwitz_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hurwitz_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ibeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ibeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ibetac" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ibetac >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_lgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_cdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_logistic_cdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_logistic_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_cdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lognormal_cdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lognormal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lpochhammer_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lpochhammer_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lpochhammer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lpochhammer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_owens_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_owens_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pochhammer_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pochhammer_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_pochhammer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_pochhammer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_cdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_cdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_psi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_psi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_qgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_qgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_radpoly" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_radpoly >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_riemann_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_shi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_shi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinc_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinhint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinhint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sinint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sin_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_sin_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "complex_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/complex_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_comp_ellint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/deathmatch_comp_ellint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_conf_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/deathmatch_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_conf_hyperg_lim" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/deathmatch_conf_hyperg_lim >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/deathmatch_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr56216_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


mpfrcalc: mpfr_gexpr.c
	$(GCC) -I. -o mpfrcalc mpfr_gexpr.c -lgmp -lmpfr -lm


test_special_function: test_special_function.cpp test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -o test_special_function test_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

diff_special_function: diff_special_function.cpp test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -o diff_special_function diff_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

testcase2: testcase2.cpp testcase2.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -o testcase2 testcase2.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

testcase: testcase.cpp testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -o testcase testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

test_limits: test_limits.cpp
	$(CXX17) -I. -o test_limits test_limits.cpp -lquadmath

test_cmath: test_cmath.cpp
	$(CXX) -o test_cmath test_cmath.cpp -lquadmath

test_airy: test_airy.cpp sf_airy.tcc
	$(CXX) -o test_airy -I$(GSL_INC_DIR) test_airy.cpp -lquadmath -Lwrappers/debug -lwrap_gsl

test_csint: test_csint.cpp csint.tcc
	$(CXX) -o test_csint test_csint.cpp -lquadmath

test_Faddeeva: $(FAD_DIR)/Faddeeva.h $(FAD_DIR)/Faddeeva.cpp
	$(CXX) -DTEST_FADDEEVA -o $(FAD_DIR)/test_Faddeeva $(FAD_DIR)/Faddeeva.cpp -lquadmath

test_fresnel: test_fresnel.cpp fresnel.tcc
	$(CXX) -o test_fresnel test_fresnel.cpp -lquadmath

test_hermite: test_hermite.cpp new_hermite.tcc
	$(CXX) -o test_hermite test_hermite.cpp -lquadmath

test_bessel: test_bessel.cpp new_bessel.tcc
	$(CXX) -o test_bessel test_bessel.cpp -lquadmath

test_nric_bessel: test_nric_bessel.cpp nric_bessel.tcc
	$(CXX) -o test_nric_bessel test_nric_bessel.cpp -lquadmath

test_anger_weber: test_anger_weber.cpp
	$(CXX17) -I. -o test_anger_weber test_anger_weber.cpp -lquadmath

test_bernoulli: test_bernoulli.cpp
	$(CXX17) -I. -o test_bernoulli test_bernoulli.cpp -lquadmath

test_bessel_asymp: test_bessel_asymp.cpp
	$(CXX17) -I. -o test_bessel_asymp test_bessel_asymp.cpp -lquadmath

test_bessel_iter: test_bessel_iter.cpp
	$(CXX17) -I. -o test_bessel_iter test_bessel_iter.cpp -lquadmath

test_beta: test_beta.cpp
	$(CXX17) -I. -o test_beta test_beta.cpp -lquadmath

test_beta_inc: test_beta_inc.cpp
	$(CXX17) -I. -o test_beta_inc test_beta_inc.cpp -lquadmath

test_binet: test_binet.cpp
	$(CXX17) -I. -o test_binet test_binet.cpp -lquadmath

test_binet_float: test_binet_float.cpp
	$(CXX17) -I. -o test_binet_float test_binet_float.cpp -lquadmath

test_bose_einstein: test_bose_einstein.cpp
	$(CXX17) -I. -o test_bose_einstein test_bose_einstein.cpp -lquadmath

test_chebyshev: test_chebyshev.cpp
	$(CXX17) -I. -o test_chebyshev test_chebyshev.cpp -lquadmath

test_chebyshev_trig: test_chebyshev_trig.cpp
	$(CXX17) -I. -o test_chebyshev_trig test_chebyshev_trig.cpp -lquadmath

test_chebyshev_trig_pi: test_chebyshev_trig_pi.cpp
	$(CXX17) -I. -o test_chebyshev_trig_pi test_chebyshev_trig_pi.cpp -lquadmath

test_clausen: test_clausen.cpp
	$(CXX17) -I. -o test_clausen test_clausen.cpp -lquadmath

test_comp_ellint_1: test_comp_ellint_1.cpp
	$(CXX17) -I. -o test_comp_ellint_1 test_comp_ellint_1.cpp -lquadmath

test_complex128: test_complex128.cpp
	$(CXX17) -I. -o test_complex128 test_complex128.cpp -lquadmath

test_complex_gamma: test_complex_gamma.cpp
	$(CXX17) -I. -o test_complex_gamma test_complex_gamma.cpp -lquadmath

test_conf_hyperg_limit: test_conf_hyperg_limit.cpp
	$(CXX17) -I. -o test_conf_hyperg_limit test_conf_hyperg_limit.cpp -lquadmath

test_const: test_const.cpp
	$(CXX17) -I. -o test_const test_const.cpp -lquadmath

test_continued_fraction: test_continued_fraction.cpp
	$(CXX17) -I. -o test_continued_fraction test_continued_fraction.cpp -lquadmath

test_cyl_hankel: test_cyl_hankel.cpp
	$(CXX17) -I. -o test_cyl_hankel test_cyl_hankel.cpp -lquadmath

test_dawson: test_dawson.cpp
	$(CXX17) -I. -o test_dawson test_dawson.cpp -lquadmath

test_debye: test_debye.cpp
	$(CXX17) -I. -o test_debye test_debye.cpp -lquadmath

test_dilog: test_dilog.cpp
	$(CXX17) -I. -o test_dilog test_dilog.cpp -lquadmath

test_expint: test_expint.cpp
	$(CXX17) -I. -o test_expint test_expint.cpp -lquadmath

test_factorial: test_factorial.cpp
	$(CXX17) -I. -o test_factorial test_factorial.cpp -lquadmath

test_fermi_dirac: test_fermi_dirac.cpp
	$(CXX17) -I. -o test_fermi_dirac test_fermi_dirac.cpp -lquadmath

test_float128: test_float128.cpp
	$(CXX17) -I. -o test_float128 test_float128.cpp -lquadmath

test_gamma: test_gamma.cpp
	$(CXX17) -I. -o test_gamma test_gamma.cpp -lquadmath

test_gamma_ratio: test_gamma_ratio.cpp
	$(CXX17) -I. -o test_gamma_ratio test_gamma_ratio.cpp -lquadmath

test_gamma_reciprocal: test_gamma_reciprocal.cpp
	$(CXX17) -I. -o test_gamma_reciprocal test_gamma_reciprocal.cpp -lquadmath

test_gegenbauer: test_gegenbauer.cpp
	$(CXX17) -I. -o test_gegenbauer test_gegenbauer.cpp -lquadmath

test_hankel: test_hankel.cpp
	$(CXX17) -I. -o test_hankel test_hankel.cpp -lquadmath

test_hankel_real_arg: test_hankel_real_arg.cpp
	$(CXX17) -I. -o test_hankel_real_arg test_hankel_real_arg.cpp -lquadmath

test_heuman_lambda: test_heuman_lambda.cpp
	$(CXX17) -I. -o test_heuman_lambda test_heuman_lambda.cpp -lquadmath

test_hurwitz_zeta: test_hurwitz_zeta.cpp
	$(CXX17) -I. -o test_hurwitz_zeta test_hurwitz_zeta.cpp -lquadmath

test_hurwitz_zeta_new: test_hurwitz_zeta_new.cpp
	$(CXX17) -I. -o test_hurwitz_zeta_new test_hurwitz_zeta_new.cpp -lquadmath

test_hyperg: test_hyperg.cpp
	$(CXX17) -I. -o test_hyperg test_hyperg.cpp -lquadmath

test_hypot: test_hypot.cpp
	$(CXX17) -I. -o test_hypot test_hypot.cpp -lquadmath

test_inv_erf: test_inv_erf.cpp
	$(CXX17) -I. -o test_inv_erf test_inv_erf.cpp -lquadmath

test_inv_ibeta: test_inv_ibeta.cpp
	$(CXX17) -I. -o test_inv_ibeta test_inv_ibeta.cpp -lquadmath

test_jacobi: test_jacobi.cpp
	$(CXX17) -I. -o test_jacobi test_jacobi.cpp -lquadmath

test_jacobi_inv: test_jacobi_inv.cpp
	$(CXX17) -I. -o test_jacobi_inv test_jacobi_inv.cpp -lquadmath

test_jacobi_zeta: test_jacobi_zeta.cpp
	$(CXX17) -I. -o test_jacobi_zeta test_jacobi_zeta.cpp -lquadmath

test_kelvin: test_kelvin.cpp
	$(CXX17) -I. -o test_kelvin test_kelvin.cpp -lquadmath

test_laguerre: test_laguerre.cpp
	$(CXX17) -I. -o test_laguerre test_laguerre.cpp -lquadmath

test_legendre: test_legendre.cpp
	$(CXX17) -I. -o test_legendre test_legendre.cpp -lquadmath

test_lentz_continued_fraction: test_lentz_continued_fraction.cpp
	$(CXX17) -I. -o test_lentz_continued_fraction test_lentz_continued_fraction.cpp -lquadmath

test_lerch: test_lerch.cpp
	$(CXX17) -I. -o test_lerch test_lerch.cpp -lquadmath

test_little_airy: test_little_airy.cpp
	$(CXX17) -I. -o test_little_airy test_little_airy.cpp -lquadmath

test_math_h: test_math_h.cpp
	$(CXX17) -I. -o test_math_h test_math_h.cpp -lquadmath

test_notsospecfun: test_notsospecfun.cpp
	$(CXX17) -I. -o test_notsospecfun test_notsospecfun.cpp -lquadmath

test_numeric_limits: test_numeric_limits.cpp
	$(CXX17) -I. -o test_numeric_limits test_numeric_limits.cpp -lquadmath

test_owens_t: test_owens_t.cpp
	$(CXX17) -I. -o test_owens_t test_owens_t.cpp -lquadmath

test_parab_cyl: test_parab_cyl.cpp
	$(CXX17) -I. -o test_parab_cyl test_parab_cyl.cpp -lquadmath

test_pochhammer: test_pochhammer.cpp
	$(CXX17) -I. -o test_pochhammer test_pochhammer.cpp -lquadmath

test_pochhammer_lower: test_pochhammer_lower.cpp
	$(CXX17) -I. -o test_pochhammer_lower test_pochhammer_lower.cpp -lquadmath

test_polylog: test_polylog.cpp
	$(CXX17) -I. -o test_polylog test_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes

test_polynomial: test_polynomial.cpp
	$(CXX17) -I. -o test_polynomial test_polynomial.cpp -lquadmath

test_pow_funs: test_pow_funs.cpp
	$(CXX17) -I. -o test_pow_funs test_pow_funs.cpp -lquadmath

test_psi: test_psi.cpp
	$(CXX17) -I. -o test_psi test_psi.cpp -lquadmath

test_rational: test_rational.cpp
	$(CXX17) -I. -o test_rational test_rational.cpp -lquadmath

test_ratpoly: test_ratpoly.cpp
	$(CXX17) -I. -o test_ratpoly test_ratpoly.cpp -lquadmath

test_recursion: test_recursion.cpp
	$(CXX17) -I. -o test_recursion test_recursion.cpp -lquadmath

test_reperiodized_hyper: test_reperiodized_hyper.cpp
	$(CXX17) -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp -lquadmath

test_reperiodized_trig: test_reperiodized_trig.cpp
	$(CXX17) -I. -o test_reperiodized_trig test_reperiodized_trig.cpp -lquadmath

test_riemann_zeta: test_riemann_zeta.cpp
	$(CXX17) -I. -o test_riemann_zeta test_riemann_zeta.cpp -lquadmath

test_root_finding: test_root_finding.cpp
	$(CXX17) -I. -o test_root_finding test_root_finding.cpp -lquadmath

test_sincos: test_sincos.cpp
	$(CXX17) -I. -o test_sincos test_sincos.cpp -lquadmath

test_sinus_cardinal: test_sinus_cardinal.cpp
	$(CXX17) -I. -o test_sinus_cardinal test_sinus_cardinal.cpp -lquadmath

test_sph_bessel: test_sph_bessel.cpp
	$(CXX17) -I. -o test_sph_bessel test_sph_bessel.cpp -lquadmath

test_sph_hankel: test_sph_hankel.cpp
	$(CXX17) -I. -o test_sph_hankel test_sph_hankel.cpp -lquadmath

test_static_polynomial: test_static_polynomial.cpp
	$(CXX17) -I. -o test_static_polynomial test_static_polynomial.cpp -lquadmath

test_steed_continued_fraction: test_steed_continued_fraction.cpp
	$(CXX17) -I. -o test_steed_continued_fraction test_steed_continued_fraction.cpp -lquadmath

test_struve: test_struve.cpp
	$(CXX17) -I. -o test_struve test_struve.cpp -lquadmath

test_struve_old: test_struve_old.cpp
	$(CXX17) -I. -o test_struve_old test_struve_old.cpp -lquadmath

test_summation: test_summation.cpp
	$(CXX17) -I. -o test_summation test_summation.cpp -lquadmath

test_theta: test_theta.cpp
	$(CXX17) -I. -o test_theta test_theta.cpp -lquadmath

test_trigamma: test_trigamma.cpp
	$(CXX17) -I. -o test_trigamma test_trigamma.cpp -lquadmath


airy_toy: airy_toy.cpp
	$(CXX17) -I. -o airy_toy airy_toy.cpp -lquadmath

hankel_toy: hankel_toy.cpp
	$(CXX17) -I. -o hankel_toy hankel_toy.cpp -lquadmath

hankel_toy128: hankel_toy128.cpp
	$(CXX17) -I. -o hankel_toy128 hankel_toy128.cpp -lquadmath

hankel_toy_new: hankel_toy_new.cpp
	$(CXX17) -I. -o hankel_toy_new hankel_toy_new.cpp -lquadmath


${CHECK_DIR}/check_airy_ai: ${CHECK_DIR}/check_airy_ai.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_ai ${CHECK_DIR}/check_airy_ai.cc -lquadmath

${CHECK_DIR}/check_airy_bi: ${CHECK_DIR}/check_airy_bi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_bi ${CHECK_DIR}/check_airy_bi.cc -lquadmath

${CHECK_DIR}/check_assoc_laguerre: ${CHECK_DIR}/check_assoc_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_laguerre ${CHECK_DIR}/check_assoc_laguerre.cc -lquadmath

${CHECK_DIR}/check_assoc_legendre: ${CHECK_DIR}/check_assoc_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_legendre ${CHECK_DIR}/check_assoc_legendre.cc -lquadmath

${CHECK_DIR}/check_bernoulli: ${CHECK_DIR}/check_bernoulli.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_bernoulli ${CHECK_DIR}/check_bernoulli.cc -lquadmath

${CHECK_DIR}/check_beta: ${CHECK_DIR}/check_beta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_beta ${CHECK_DIR}/check_beta.cc -lquadmath

${CHECK_DIR}/check_bincoef: ${CHECK_DIR}/check_bincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_bincoef ${CHECK_DIR}/check_bincoef.cc -lquadmath

${CHECK_DIR}/check_chebyshev_t: ${CHECK_DIR}/check_chebyshev_t.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_t ${CHECK_DIR}/check_chebyshev_t.cc -lquadmath

${CHECK_DIR}/check_chebyshev_u: ${CHECK_DIR}/check_chebyshev_u.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_u ${CHECK_DIR}/check_chebyshev_u.cc -lquadmath

${CHECK_DIR}/check_chebyshev_v: ${CHECK_DIR}/check_chebyshev_v.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_v ${CHECK_DIR}/check_chebyshev_v.cc -lquadmath

${CHECK_DIR}/check_chebyshev_w: ${CHECK_DIR}/check_chebyshev_w.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_w ${CHECK_DIR}/check_chebyshev_w.cc -lquadmath

${CHECK_DIR}/check_chi: ${CHECK_DIR}/check_chi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chi ${CHECK_DIR}/check_chi.cc -lquadmath

${CHECK_DIR}/check_clausen_c: ${CHECK_DIR}/check_clausen_c.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_clausen_c ${CHECK_DIR}/check_clausen_c.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_1: ${CHECK_DIR}/check_comp_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_1 ${CHECK_DIR}/check_comp_ellint_1.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_2: ${CHECK_DIR}/check_comp_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_2 ${CHECK_DIR}/check_comp_ellint_2.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_3: ${CHECK_DIR}/check_comp_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_3 ${CHECK_DIR}/check_comp_ellint_3.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_d: ${CHECK_DIR}/check_comp_ellint_d.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_d ${CHECK_DIR}/check_comp_ellint_d.cc -lquadmath

${CHECK_DIR}/check_conf_hyperg: ${CHECK_DIR}/check_conf_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg ${CHECK_DIR}/check_conf_hyperg.cc -lquadmath

${CHECK_DIR}/check_conf_hyperg_lim: ${CHECK_DIR}/check_conf_hyperg_lim.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg_lim ${CHECK_DIR}/check_conf_hyperg_lim.cc -lquadmath

${CHECK_DIR}/check_coshint: ${CHECK_DIR}/check_coshint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_coshint ${CHECK_DIR}/check_coshint.cc -lquadmath

${CHECK_DIR}/check_cosint: ${CHECK_DIR}/check_cosint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cosint ${CHECK_DIR}/check_cosint.cc -lquadmath

${CHECK_DIR}/check_cos_pi: ${CHECK_DIR}/check_cos_pi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cos_pi ${CHECK_DIR}/check_cos_pi.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_i: ${CHECK_DIR}/check_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_i ${CHECK_DIR}/check_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_j: ${CHECK_DIR}/check_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_j ${CHECK_DIR}/check_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_k: ${CHECK_DIR}/check_cyl_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_k ${CHECK_DIR}/check_cyl_bessel_k.cc -lquadmath

${CHECK_DIR}/check_cyl_hankel_1: ${CHECK_DIR}/check_cyl_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_1 ${CHECK_DIR}/check_cyl_hankel_1.cc -lquadmath

${CHECK_DIR}/check_cyl_hankel_2: ${CHECK_DIR}/check_cyl_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_2 ${CHECK_DIR}/check_cyl_hankel_2.cc -lquadmath

${CHECK_DIR}/check_cyl_neumann: ${CHECK_DIR}/check_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_neumann ${CHECK_DIR}/check_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_dawson: ${CHECK_DIR}/check_dawson.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dawson ${CHECK_DIR}/check_dawson.cc -lquadmath -lquadmath

${CHECK_DIR}/check_dilog: ${CHECK_DIR}/check_dilog.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dilog ${CHECK_DIR}/check_dilog.cc -lquadmath -lquadmath

${CHECK_DIR}/check_dirichlet_beta: ${CHECK_DIR}/check_dirichlet_beta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_beta ${CHECK_DIR}/check_dirichlet_beta.cc -lquadmath

${CHECK_DIR}/check_dirichlet_eta: ${CHECK_DIR}/check_dirichlet_eta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_eta ${CHECK_DIR}/check_dirichlet_eta.cc -lquadmath

${CHECK_DIR}/check_dirichlet_lambda: ${CHECK_DIR}/check_dirichlet_lambda.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_lambda ${CHECK_DIR}/check_dirichlet_lambda.cc -lquadmath

${CHECK_DIR}/check_double_factorial: ${CHECK_DIR}/check_double_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_double_factorial ${CHECK_DIR}/check_double_factorial.cc -lquadmath

${CHECK_DIR}/check_ellint_1: ${CHECK_DIR}/check_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_1 ${CHECK_DIR}/check_ellint_1.cc -lquadmath

${CHECK_DIR}/check_ellint_2: ${CHECK_DIR}/check_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_2 ${CHECK_DIR}/check_ellint_2.cc -lquadmath

${CHECK_DIR}/check_ellint_3: ${CHECK_DIR}/check_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_3 ${CHECK_DIR}/check_ellint_3.cc -lquadmath

${CHECK_DIR}/check_ellint_d: ${CHECK_DIR}/check_ellint_d.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_d ${CHECK_DIR}/check_ellint_d.cc -lquadmath

${CHECK_DIR}/check_ellint_rc: ${CHECK_DIR}/check_ellint_rc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rc ${CHECK_DIR}/check_ellint_rc.cc -lquadmath

${CHECK_DIR}/check_ellint_rd: ${CHECK_DIR}/check_ellint_rd.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rd ${CHECK_DIR}/check_ellint_rd.cc -lquadmath

${CHECK_DIR}/check_ellint_rf: ${CHECK_DIR}/check_ellint_rf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rf ${CHECK_DIR}/check_ellint_rf.cc -lquadmath

${CHECK_DIR}/check_ellint_rg: ${CHECK_DIR}/check_ellint_rg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rg ${CHECK_DIR}/check_ellint_rg.cc -lquadmath

${CHECK_DIR}/check_ellint_rj: ${CHECK_DIR}/check_ellint_rj.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rj ${CHECK_DIR}/check_ellint_rj.cc -lquadmath

${CHECK_DIR}/check_ellnome: ${CHECK_DIR}/check_ellnome.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellnome ${CHECK_DIR}/check_ellnome.cc -lquadmath

${CHECK_DIR}/check_expint: ${CHECK_DIR}/check_expint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint ${CHECK_DIR}/check_expint.cc -lquadmath

${CHECK_DIR}/check_expint_en: ${CHECK_DIR}/check_expint_en.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint_en ${CHECK_DIR}/check_expint_en.cc -lquadmath

${CHECK_DIR}/check_factorial: ${CHECK_DIR}/check_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_factorial ${CHECK_DIR}/check_factorial.cc -lquadmath

${CHECK_DIR}/check_fresnel_c: ${CHECK_DIR}/check_fresnel_c.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_c ${CHECK_DIR}/check_fresnel_c.cc -lquadmath

${CHECK_DIR}/check_fresnel_s: ${CHECK_DIR}/check_fresnel_s.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_s ${CHECK_DIR}/check_fresnel_s.cc -lquadmath

${CHECK_DIR}/check_gamma_reciprocal: ${CHECK_DIR}/check_gamma_reciprocal.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_reciprocal ${CHECK_DIR}/check_gamma_reciprocal.cc -lquadmath

${CHECK_DIR}/check_gegenbauer: ${CHECK_DIR}/check_gegenbauer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gegenbauer ${CHECK_DIR}/check_gegenbauer.cc -lquadmath

${CHECK_DIR}/check_hermite: ${CHECK_DIR}/check_hermite.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hermite ${CHECK_DIR}/check_hermite.cc -lquadmath

${CHECK_DIR}/check_heuman_lambda: ${CHECK_DIR}/check_heuman_lambda.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_heuman_lambda ${CHECK_DIR}/check_heuman_lambda.cc -lquadmath

${CHECK_DIR}/check_hurwitz_zeta: ${CHECK_DIR}/check_hurwitz_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hurwitz_zeta ${CHECK_DIR}/check_hurwitz_zeta.cc -lquadmath

${CHECK_DIR}/check_hyperg: ${CHECK_DIR}/check_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hyperg ${CHECK_DIR}/check_hyperg.cc -lquadmath

${CHECK_DIR}/check_ibeta: ${CHECK_DIR}/check_ibeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ibeta ${CHECK_DIR}/check_ibeta.cc -lquadmath

${CHECK_DIR}/check_ibetac: ${CHECK_DIR}/check_ibetac.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ibetac ${CHECK_DIR}/check_ibetac.cc -lquadmath

${CHECK_DIR}/check_jacobi: ${CHECK_DIR}/check_jacobi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi ${CHECK_DIR}/check_jacobi.cc -lquadmath

${CHECK_DIR}/check_jacobi_cn: ${CHECK_DIR}/check_jacobi_cn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_cn ${CHECK_DIR}/check_jacobi_cn.cc -lquadmath

${CHECK_DIR}/check_jacobi_dn: ${CHECK_DIR}/check_jacobi_dn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_dn ${CHECK_DIR}/check_jacobi_dn.cc -lquadmath

${CHECK_DIR}/check_jacobi_sn: ${CHECK_DIR}/check_jacobi_sn.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_sn ${CHECK_DIR}/check_jacobi_sn.cc -lquadmath

${CHECK_DIR}/check_jacobi_zeta: ${CHECK_DIR}/check_jacobi_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_zeta ${CHECK_DIR}/check_jacobi_zeta.cc -lquadmath

${CHECK_DIR}/check_laguerre: ${CHECK_DIR}/check_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_laguerre ${CHECK_DIR}/check_laguerre.cc -lquadmath

${CHECK_DIR}/check_lbincoef: ${CHECK_DIR}/check_lbincoef.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lbincoef ${CHECK_DIR}/check_lbincoef.cc -lquadmath

${CHECK_DIR}/check_ldouble_factorial: ${CHECK_DIR}/check_ldouble_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ldouble_factorial ${CHECK_DIR}/check_ldouble_factorial.cc -lquadmath

${CHECK_DIR}/check_legendre: ${CHECK_DIR}/check_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre ${CHECK_DIR}/check_legendre.cc -lquadmath

${CHECK_DIR}/check_legendre_q: ${CHECK_DIR}/check_legendre_q.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre_q ${CHECK_DIR}/check_legendre_q.cc -lquadmath

${CHECK_DIR}/check_lfactorial: ${CHECK_DIR}/check_lfactorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfactorial ${CHECK_DIR}/check_lfactorial.cc -lquadmath

${CHECK_DIR}/check_lgamma: ${CHECK_DIR}/check_lgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lgamma ${CHECK_DIR}/check_lgamma.cc -lquadmath

${CHECK_DIR}/check_logistic_cdf: ${CHECK_DIR}/check_logistic_cdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_cdf ${CHECK_DIR}/check_logistic_cdf.cc -lquadmath

${CHECK_DIR}/check_logistic_pdf: ${CHECK_DIR}/check_logistic_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_pdf ${CHECK_DIR}/check_logistic_pdf.cc -lquadmath

${CHECK_DIR}/check_lognormal_cdf: ${CHECK_DIR}/check_lognormal_cdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_cdf ${CHECK_DIR}/check_lognormal_cdf.cc -lquadmath

${CHECK_DIR}/check_lognormal_pdf: ${CHECK_DIR}/check_lognormal_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_pdf ${CHECK_DIR}/check_lognormal_pdf.cc -lquadmath

${CHECK_DIR}/check_lpochhammer_lower: ${CHECK_DIR}/check_lpochhammer_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lpochhammer_lower ${CHECK_DIR}/check_lpochhammer_lower.cc -lquadmath

${CHECK_DIR}/check_lpochhammer: ${CHECK_DIR}/check_lpochhammer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lpochhammer ${CHECK_DIR}/check_lpochhammer.cc -lquadmath

${CHECK_DIR}/check_normal_cdf: ${CHECK_DIR}/check_normal_cdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_cdf ${CHECK_DIR}/check_normal_cdf.cc -lquadmath

${CHECK_DIR}/check_normal_pdf: ${CHECK_DIR}/check_normal_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_pdf ${CHECK_DIR}/check_normal_pdf.cc -lquadmath

${CHECK_DIR}/check_owens_t: ${CHECK_DIR}/check_owens_t.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_owens_t ${CHECK_DIR}/check_owens_t.cc -lquadmath

${CHECK_DIR}/check_pgamma: ${CHECK_DIR}/check_pgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pgamma ${CHECK_DIR}/check_pgamma.cc -lquadmath

${CHECK_DIR}/check_pochhammer_lower: ${CHECK_DIR}/check_pochhammer_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pochhammer_lower ${CHECK_DIR}/check_pochhammer_lower.cc -lquadmath

${CHECK_DIR}/check_pochhammer: ${CHECK_DIR}/check_pochhammer.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_pochhammer ${CHECK_DIR}/check_pochhammer.cc -lquadmath

${CHECK_DIR}/check_psi: ${CHECK_DIR}/check_psi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_psi ${CHECK_DIR}/check_psi.cc -lquadmath

${CHECK_DIR}/check_qgamma: ${CHECK_DIR}/check_qgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_qgamma ${CHECK_DIR}/check_qgamma.cc -lquadmath

${CHECK_DIR}/check_radpoly: ${CHECK_DIR}/check_radpoly.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_radpoly ${CHECK_DIR}/check_radpoly.cc -lquadmath

${CHECK_DIR}/check_riemann_zeta: ${CHECK_DIR}/check_riemann_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_riemann_zeta ${CHECK_DIR}/check_riemann_zeta.cc -lquadmath

${CHECK_DIR}/check_shi: ${CHECK_DIR}/check_shi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_shi ${CHECK_DIR}/check_shi.cc -lquadmath

${CHECK_DIR}/check_sinc: ${CHECK_DIR}/check_sinc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc ${CHECK_DIR}/check_sinc.cc -lquadmath

${CHECK_DIR}/check_sinc_pi: ${CHECK_DIR}/check_sinc_pi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc_pi ${CHECK_DIR}/check_sinc_pi.cc -lquadmath

${CHECK_DIR}/check_sinhint: ${CHECK_DIR}/check_sinhint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinhint ${CHECK_DIR}/check_sinhint.cc -lquadmath

${CHECK_DIR}/check_sinint: ${CHECK_DIR}/check_sinint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinint ${CHECK_DIR}/check_sinint.cc -lquadmath

${CHECK_DIR}/check_sin_pi: ${CHECK_DIR}/check_sin_pi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sin_pi ${CHECK_DIR}/check_sin_pi.cc -lquadmath

${CHECK_DIR}/check_sph_bessel: ${CHECK_DIR}/check_sph_bessel.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel ${CHECK_DIR}/check_sph_bessel.cc -lquadmath

${CHECK_DIR}/check_sph_bessel_i: ${CHECK_DIR}/check_sph_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_i ${CHECK_DIR}/check_sph_bessel_i.cc -lquadmath

${CHECK_DIR}/check_sph_bessel_k: ${CHECK_DIR}/check_sph_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_k ${CHECK_DIR}/check_sph_bessel_k.cc -lquadmath

${CHECK_DIR}/check_sph_hankel_1: ${CHECK_DIR}/check_sph_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_1 ${CHECK_DIR}/check_sph_hankel_1.cc -lquadmath

${CHECK_DIR}/check_sph_hankel_2: ${CHECK_DIR}/check_sph_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_2 ${CHECK_DIR}/check_sph_hankel_2.cc -lquadmath

${CHECK_DIR}/check_sph_harmonic: ${CHECK_DIR}/check_sph_harmonic.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_harmonic ${CHECK_DIR}/check_sph_harmonic.cc -lquadmath

${CHECK_DIR}/check_sph_legendre: ${CHECK_DIR}/check_sph_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_legendre ${CHECK_DIR}/check_sph_legendre.cc -lquadmath

${CHECK_DIR}/check_sph_neumann: ${CHECK_DIR}/check_sph_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_neumann ${CHECK_DIR}/check_sph_neumann.cc -lquadmath

${CHECK_DIR}/check_tgamma_lower: ${CHECK_DIR}/check_tgamma_lower.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma_lower ${CHECK_DIR}/check_tgamma_lower.cc -lquadmath

${CHECK_DIR}/check_tgamma: ${CHECK_DIR}/check_tgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma ${CHECK_DIR}/check_tgamma.cc -lquadmath

${CHECK_DIR}/check_theta_1: ${CHECK_DIR}/check_theta_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_1 ${CHECK_DIR}/check_theta_1.cc -lquadmath

${CHECK_DIR}/check_theta_2: ${CHECK_DIR}/check_theta_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_2 ${CHECK_DIR}/check_theta_2.cc -lquadmath

${CHECK_DIR}/check_theta_3: ${CHECK_DIR}/check_theta_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_3 ${CHECK_DIR}/check_theta_3.cc -lquadmath

${CHECK_DIR}/check_theta_4: ${CHECK_DIR}/check_theta_4.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_4 ${CHECK_DIR}/check_theta_4.cc -lquadmath

${CHECK_DIR}/check_zernike: ${CHECK_DIR}/check_zernike.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_zernike ${CHECK_DIR}/check_zernike.cc -lquadmath

${CHECK_DIR}/complex_ellint_rc: ${CHECK_DIR}/complex_ellint_rc.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rc ${CHECK_DIR}/complex_ellint_rc.cc -lquadmath

${CHECK_DIR}/complex_ellint_rd: ${CHECK_DIR}/complex_ellint_rd.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rd ${CHECK_DIR}/complex_ellint_rd.cc -lquadmath

${CHECK_DIR}/complex_ellint_rf: ${CHECK_DIR}/complex_ellint_rf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rf ${CHECK_DIR}/complex_ellint_rf.cc -lquadmath

${CHECK_DIR}/complex_ellint_rg: ${CHECK_DIR}/complex_ellint_rg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rg ${CHECK_DIR}/complex_ellint_rg.cc -lquadmath

${CHECK_DIR}/complex_ellint_rj: ${CHECK_DIR}/complex_ellint_rj.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rj ${CHECK_DIR}/complex_ellint_rj.cc -lquadmath

${CHECK_DIR}/complex_airy_ai: ${CHECK_DIR}/complex_airy_ai.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_airy_ai ${CHECK_DIR}/complex_airy_ai.cc -lquadmath

${CHECK_DIR}/complex_airy_bi: ${CHECK_DIR}/complex_airy_bi.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_airy_bi ${CHECK_DIR}/complex_airy_bi.cc -lquadmath

${CHECK_DIR}/deathmatch_comp_ellint: ${CHECK_DIR}/deathmatch_comp_ellint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_comp_ellint ${CHECK_DIR}/deathmatch_comp_ellint.cc -lquadmath

${CHECK_DIR}/deathmatch_conf_hyperg: ${CHECK_DIR}/deathmatch_conf_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_conf_hyperg ${CHECK_DIR}/deathmatch_conf_hyperg.cc -lquadmath

${CHECK_DIR}/deathmatch_conf_hyperg_lim: ${CHECK_DIR}/deathmatch_conf_hyperg_lim.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_conf_hyperg_lim ${CHECK_DIR}/deathmatch_conf_hyperg_lim.cc -lquadmath

${CHECK_DIR}/deathmatch_hyperg: ${CHECK_DIR}/deathmatch_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_hyperg ${CHECK_DIR}/deathmatch_hyperg.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_hankel_1: ${CHECK_DIR}/pr56216_cyl_hankel_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_1 ${CHECK_DIR}/pr56216_cyl_hankel_1.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_hankel_2: ${CHECK_DIR}/pr56216_cyl_hankel_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_2 ${CHECK_DIR}/pr56216_cyl_hankel_2.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_bessel_i: ${CHECK_DIR}/pr56216_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_bessel_i ${CHECK_DIR}/pr56216_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/origin_cyl_bessel_j: ${CHECK_DIR}/origin_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_bessel_j ${CHECK_DIR}/origin_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/origin_cyl_neumann: ${CHECK_DIR}/origin_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_neumann ${CHECK_DIR}/origin_cyl_neumann.cc -lquadmath


tarball:
	mkdir tr29124
	cp Makefile cmath *.* tr29124
	tar -cvf tr29124.tar tr29124
	bzip2 -f tr29124.tar
	md5sum tr29124.tar.bz2 > tr29124.tar.bz2.md5
	rm -rf tr29124

clean:
	rm -f $(BINS)
	rm -f $(CHECKS)
	rm -f tr29124.tar.bz2 tr29124.tar.bz2.md5
	rm -f *.stackdump

testclean:
	rm -f test/gsl_*_[fdl].txt
	rm -f test/std_*_[fdl].txt
	rm -f test/tr1_*_[fdl].txt

diffclean:
	rm -f diff/diff_*_[fdl].txt
