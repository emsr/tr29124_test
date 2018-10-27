
# -Wconversion

# Don't look for this gcc anymore.
#SUFFIX = _tr29124
# This will start with $(HOME)/bin...
CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
ifeq ("$(wildcard $(CXX_INST_DIR))","")
  SUFFIX = 
  CXX_INST_DIR = $(HOME)/bin
  ifeq ("$(wildcard $(CXX_INST_DIR))","")
    ifneq ($(wildcard "/mingw64"),"")
      CXX_INST_DIR = /mingw64
    else
      CXX_INST_DIR = /usr
    endif
  endif
endif

#  -fsanitize=address
GFORTRAN = $(CXX_INST_DIR)/bin/gfortran -g -Wall -Wextra -Wno-compare-reals
GCC = $(CXX_INST_DIR)/bin/gcc -g -Wall -Wextra
CXX = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi
CXX_INC_DIR = $(CXX_INST_DIR)/include/c++/7.0.0/bits
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = .

INC_DIR = bits

FAD_DIR = Faddeeva

MATH_DIR = ./libstdc++_support

CHECK_DIR = $(MATH_DIR)/check

OBJ_DIR = $(MATH_DIR)/obj

LIB_DIR = $(MATH_DIR)


#       test_special_function \
#       testcase2 \
#

BINS = libstdc++_support/testcase \
       libstdc++_support/testcase_tr1 \
       mpfrcalc \
       diff_special_function \
       test_special_function \
       airy_toy \
       hankel_toy \
       hankel_toy128 \
       hankel_toy_new \
       test_airy_roots \
       test_anger_weber \
       test_appell_f1 \
       test_arith_geom_mean \
       test_bernoulli \
       test_bessel \
       test_bessel_asymp \
       test_bessel_iter \
       test_beta \
       test_beta_inc \
       test_binet \
       test_binet_float \
       test_bose_einstein \
       test_charlier \
       test_chebyshev \
       test_chebyshev_trig \
       test_chebyshev_trig_pi \
       test_clausen \
       test_cmath \
       test_comp_ellint_1 \
       test_complex128 \
       test_complex_gamma \
       test_conf_hyperg \
       test_conf_hyperg_limit \
       test_const \
       test_continued_fraction \
       test_continuous_dual_hahn \
       test_continuous_hahn \
       test_cordic \
       test_coulomb \
       test_csint \
       test_cyl_hankel \
       test_dawson \
       test_debye \
       test_digamma \
       test_dilog \
       test_dirichlet_eta \
       test_dual_hahn \
       test_erfc \
       test_experfc \
       test_expint \
       test_factorial \
       test_faddeeva \
       test_falling_factorial \
       test_fermi_dirac \
       test_fibonacci \
       test_float128 \
       test_fresnel \
       test_gamma \
       test_gamma_ratio \
       test_gamma_reciprocal \
       test_gegenbauer \
       test_gudermannian \
       test_hahn \
       test_hankel \
       test_hankel_real_arg \
       test_hermite \
       test_heuman_lambda \
       test_hurwitz_zeta \
       test_hurwitz_zeta_new \
       test_hydrogen \
       test_hyperg \
       test_hypot \
       test_inv_erf \
       test_inv_gamma \
       test_inv_ibeta \
       test_jacobi \
       test_jacobi_ellint \
       test_jacobi_inv \
       test_jacobi_theta \
       test_jacobi_zeta \
       test_kelvin \
       test_krawtchouk \
       test_laguerre \
       test_lambert_w \
       test_large_order_bessel \
       test_legendre \
       test_legendre_ellint \
       test_lentz_continued_fraction \
       test_lerch \
       test_limits \
       test_little_airy \
       test_lobatto \
       test_logsumexp \
       test_lommel \
       test_marcum_q \
       test_math_h \
       test_maxint \
       test_meixner \
       test_meixner_pollaczek \
       test_mittag_leffler \
       test_mod2pi \
       test_mpreal \
       test_notsospecfun \
       test_nric_bessel \
       test_numeric_limits \
       test_owens_t \
       test_parab_cyl \
       test_polygamma \
       test_polylog \
       test_power_mean \
       test_power_norm \
       test_racah \
       test_rational \
       test_recursion \
       test_reperiodized_hyper \
       test_reperiodized_trig \
       test_riemann_zeta \
       test_rising_factorial \
       test_root_finding \
       test_sincos \
       test_sinus_cardinal \
       test_sph_bessel \
       test_sph_hankel \
       test_steed_continued_fraction \
       test_struve \
       test_struve_old \
       test_summation \
       test_theta \
       test_tr1_cmath \
       test_tricomi_u \
       test_trig \
       test_weierstrass_ellint \
       test_wilson \
       test_wright_omega \
       test_zeta_trig \
       run_coulfg \
       RUN_COULFG


CHECKS = ${CHECK_DIR}/check_airy_ai \
	 ${CHECK_DIR}/check_airy_bi \
	 ${CHECK_DIR}/check_assoc_laguerre \
	 ${CHECK_DIR}/check_assoc_legendre \
	 ${CHECK_DIR}/check_beta \
	 ${CHECK_DIR}/check_bernoulli \
	 ${CHECK_DIR}/check_binomial \
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
	 ${CHECK_DIR}/check_digamma \
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
	 ${CHECK_DIR}/check_euler \
	 ${CHECK_DIR}/check_eulerian_1 \
	 ${CHECK_DIR}/check_eulerian_2 \
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
	 ${CHECK_DIR}/check_lbinomial \
	 ${CHECK_DIR}/check_ldouble_factorial \
	 ${CHECK_DIR}/check_legendre \
	 ${CHECK_DIR}/check_legendre \
	 ${CHECK_DIR}/check_legendre_q \
	 ${CHECK_DIR}/check_lfactorial \
	 ${CHECK_DIR}/check_lgamma \
	 ${CHECK_DIR}/check_logistic_p \
	 ${CHECK_DIR}/check_logistic_pdf \
	 ${CHECK_DIR}/check_lfalling_factorial \
	 ${CHECK_DIR}/check_lrising_factorial \
	 ${CHECK_DIR}/check_lognormal_p \
	 ${CHECK_DIR}/check_lognormal_pdf \
	 ${CHECK_DIR}/check_normal_p \
	 ${CHECK_DIR}/check_normal_pdf \
	 ${CHECK_DIR}/check_owens_t \
	 ${CHECK_DIR}/check_gamma_p \
	 ${CHECK_DIR}/check_falling_factorial \
	 ${CHECK_DIR}/check_rising_factorial \
	 ${CHECK_DIR}/check_gamma_q \
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
	 ${CHECK_DIR}/check_stirling_1 \
	 ${CHECK_DIR}/check_stirling_2 \
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
	 ${CHECK_DIR}/check_clausen_cl \
	 ${CHECK_DIR}/pr56216_cyl_hankel_1 \
	 ${CHECK_DIR}/pr56216_cyl_hankel_2 \
	 ${CHECK_DIR}/pr56216_cyl_bessel_i \
	 ${CHECK_DIR}/pr68397 \
	 ${CHECK_DIR}/origin_cyl_bessel_j \
	 ${CHECK_DIR}/origin_cyl_neumann

TR1_CHECKS =  \
	${CHECK_DIR}/check_tr1_assoc_laguerre \
	${CHECK_DIR}/check_tr1_assoc_legendre \
	${CHECK_DIR}/check_tr1_beta \
	${CHECK_DIR}/check_tr1_comp_ellint_1 \
	${CHECK_DIR}/check_tr1_comp_ellint_2 \
	${CHECK_DIR}/check_tr1_comp_ellint_3 \
	${CHECK_DIR}/check_tr1_conf_hyperg \
	${CHECK_DIR}/check_tr1_cyl_bessel_i \
	${CHECK_DIR}/check_tr1_cyl_bessel_j \
	${CHECK_DIR}/check_tr1_cyl_bessel_k \
	${CHECK_DIR}/check_tr1_cyl_neumann \
	${CHECK_DIR}/check_tr1_ellint_1 \
	${CHECK_DIR}/check_tr1_ellint_2 \
	${CHECK_DIR}/check_tr1_ellint_3 \
	${CHECK_DIR}/check_tr1_expint \
	${CHECK_DIR}/check_tr1_hermite \
	${CHECK_DIR}/check_tr1_hyperg \
	${CHECK_DIR}/check_tr1_laguerre \
	${CHECK_DIR}/check_tr1_legendre \
	${CHECK_DIR}/check_tr1_riemann_zeta \
	${CHECK_DIR}/check_tr1_sph_bessel \
	${CHECK_DIR}/check_tr1_sph_legendre \
	${CHECK_DIR}/check_tr1_sph_neumann


all: $(BINS)


WRAP_DIR = wrappers
WRAP_DEBUG_DIR = $(WRAP_DIR)/debug
WRAP_RELEASE_DIR = $(WRAP_DIR)/release

WRAPPER_INCS = \
  $(WRAP_DIR)/wrap_boost.h \
  $(WRAP_DIR)/wrap_burkhardt.h \
  $(WRAP_DIR)/wrap_cephes.h \
  $(WRAP_DIR)/wrap_faddeeva.h \
  $(WRAP_DIR)/wrap_gsl.h \
  $(WRAP_DIR)/wrap_lambert.h \
  $(WRAP_DIR)/wrap_lerch.h

WRAPPER_SRCS = \
  $(WRAP_DIR)/wrap_boost.cpp \
  $(WRAP_DIR)/wrap_burkhardt.cpp \
  $(WRAP_DIR)/wrap_cephes.cpp \
  $(WRAP_DIR)/wrap_faddeeva.cpp \
  $(WRAP_DIR)/wrap_gsl.cpp \
  $(WRAP_DIR)/wrap_lambert.cpp \
  $(WRAP_DIR)/wrap_lerch.cpp

WRAP_DEBUG_LIBS = \
  $(WRAP_DEBUG_DIR)/libwrap_boost.so \
  $(WRAP_DEBUG_DIR)/libwrap_burkhardt.so \
  $(WRAP_DEBUG_DIR)/libwrap_cephes.so \
  $(WRAP_DEBUG_DIR)/libwrap_faddeeva.so \
  $(WRAP_DEBUG_DIR)/libwrap_gsl.so \
  $(WRAP_DEBUG_DIR)/libwrap_lambert.so \
  $(WRAP_DEBUG_DIR)/libwrap_lerchphi.so

WRAP_RELEASE_LIBS = \
  $(WRAP_RELEASE_DIR)/libwrap_boost.so \
  $(WRAP_RELEASE_DIR)/libwrap_burkhardt.so \
  $(WRAP_RELEASE_DIR)/libwrap_cephes.so \
  $(WRAP_RELEASE_DIR)/libwrap_faddeeva.so \
  $(WRAP_RELEASE_DIR)/libwrap_gsl.so \
  $(WRAP_RELEASE_DIR)/libwrap_lambert.so \
  $(WRAP_RELEASE_DIR)/libwrap_lerchphi.so

wrappers_debug: $(WRAP_DEBUG_DIR) $(WRAP_DEBUG_LIBS)
	touch wrappers_debug

$(WRAP_DEBUG_LIBS): $(WRAP_DEBUG_DIR)/Makefile $(WRAPPER_INCS) $(WRAPPER_SRCS)
	$(MAKE) -C $(WRAP_DEBUG_DIR)

$(WRAP_DEBUG_DIR)/Makefile: $(WRAP_DIR)/CMakeLists.txt
	(cd $(WRAP_DEBUG_DIR); rm -rf *; cmake .. -DCMAKE_BUILD_TYPE=DEBUG -G"Unix Makefiles")

$(WRAP_DEBUG_DIR):
	if test ! -d $(WRAP_DEBUG_DIR); then \
	  mkdir -p $(WRAP_DEBUG_DIR); \
	fi

wrappers_release: $(WRAP_RELEASE_DIR) $(WRAP_RELEASE_LIBS)
	touch wrappers_release

$(WRAP_RELEASE_LIBS): $(WRAP_RELEASE_DIR)/Makefile $(WRAPPER_INCS) $(WRAPPER_SRCS)
	$(MAKE) -C $(WRAP_RELEASE_DIR)

$(WRAP_RELEASE_DIR)/Makefile: $(WRAP_DIR)/CMakeLists.txt
	(cd $(WRAP_RELEASE_DIR); rm -rf *; cmake .. -DCMAKE_BUILD_TYPE=RELEASE -G"Unix Makefiles")

$(WRAP_RELEASE_DIR):
	if test ! -d $(WRAP_RELEASE_DIR); then \
	  mkdir -p $(WRAP_RELEASE_DIR); \
	fi

LambertW:
	ifeq ("$(wildcard LambertW)",""); then \
	  git clone https://github.com/emsr/LambertW.git; \
	endif

mpreal:
	cd $(HOME) && hg clone https://bitbucket.org/advanpix/mpreal


docs: bits/*
	rm -rf html/*
	rm -rf latex/*
	doxygen
	cd latex && make

testcases2: testcase2
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./testcase2

testcases: testcase
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./testcase

testcases_tr1: testcase_tr1
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./testcase_tr1

diffs: diff_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./diff_special_function > diff_special_function.txt

tests: test_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_special_function > test_special_function.txt

# This will always build the executables!
test: $(BINS)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./airy_toy > airy_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy > hankel_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy128 > hankel_toy128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./hankel_toy_new > hankel_toy_new.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_airy_roots > test_airy_roots.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_anger_weber > test_anger_weber.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_appell_f1 > test_appell_f1.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_bernoulli > test_bernoulli.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bessel > test_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bessel_asymp > test_bessel_asymp.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bessel_iter > test_bessel_iter.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_beta > test_beta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_beta_inc > test_beta_inc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_binet > test_binet.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_binet_float > test_binet_float.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_bose_einstein > test_bose_einstein.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_charlier > test_charlier.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_chebyshev > test_chebyshev.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_chebyshev_trig > test_chebyshev_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_chebyshev_trig_pi > test_chebyshev_trig_pi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_clausen > test_clausen.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_cmath > test_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_comp_ellint_1 > test_comp_ellint_1.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_complex128 > test_complex128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_complex_gamma > test_complex_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_conf_hyperg > test_conf_hyperg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_conf_hyperg_limit > test_conf_hyperg_limit.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_const > test_const.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_continued_fraction > test_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_continuous_dual_hahn > test_continuous_dual_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_continuous_hahn > test_continuous_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_cordic > test_cordic.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_coulomb > test_coulomb.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_csint > test_csint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_cyl_hankel > test_cyl_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_dawson > test_dawson.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_debye > test_debye.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_digamma > test_digamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_dilog > test_dilog.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_dirichlet_eta > test_dirichlet_eta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_dual_hahn > test_dual_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_erfc > test_erfc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_experfc > test_experfc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_expint > test_expint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_factorial > test_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_faddeeva > test_faddeeva.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_falling_factorial > test_falling_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_fermi_dirac > test_fermi_dirac.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_fibonacci > test_fibonacci.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_float128 > test_float128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_fresnel > test_fresnel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_gamma > test_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_gamma_ratio > test_gamma_ratio.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_gamma_reciprocal > test_gamma_reciprocal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_gegenbauer > test_gegenbauer.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_gudermannian > test_gudermannian.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hahn > test_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hankel > test_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_hankel_real_arg > test_hankel_real_arg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hermite > test_hermite.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_heuman_lambda > test_heuman_lambda.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hurwitz_zeta > test_hurwitz_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hurwitz_zeta_new > test_hurwitz_zeta_new.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_hydrogen > test_hydrogen.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_hyperg > test_hyperg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_hypot > test_hypot.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_inv_gamma > test_inv_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_inv_ibeta > test_inv_ibeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_jacobi > test_jacobi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_jacobi_ellint > test_jacobi_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_jacobi_inv > test_jacobi_inv.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_jacobi_theta > test_jacobi_theta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_jacobi_zeta > test_jacobi_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_kelvin > test_kelvin.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_krawtchouk > test_krawtchouk.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_laguerre > test_laguerre.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_lambert_w > test_lambert_w.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_legendre > test_legendre.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_legendre_ellint > test_legendre_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_lentz_continued_fraction > test_lentz_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_lerch > test_lerch.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_limits > test_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_little_airy > test_little_airy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_lobatto > test_lobatto.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_logsumexp > test_logsumexp.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_lommel > test_lommel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_marcum_q > test_marcum_q.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_math_h > test_math_h.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_maxint > test_maxint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_meixner > test_meixner.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_meixner_pollaczek > test_meixner_pollaczek.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_mittag_leffler > test_mittag_leffler.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_mod2pi > test_mod2pi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_mpreal > test_mpreal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_notsospecfun > test_notsospecfun.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_nric_bessel > test_nric_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_numeric_limits > test_numeric_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_owens_t > test_owens_t.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_parab_cyl > test_parab_cyl.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_polygamma > test_polygamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_polylog > test_polylog.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_power_mean > test_power_mean.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_power_norm > test_power_norm.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_racah > test_racah.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_rational > test_rational.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_recursion > test_recursion.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_reperiodized_hyper > test_reperiodized_hyper.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_reperiodized_trig > test_reperiodized_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_riemann_zeta > test_riemann_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_rising_factorial > test_rising_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_root_finding > test_root_finding.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_sincos > test_sincos.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_sinus_cardinal > test_sinus_cardinal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_sph_bessel > test_sph_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_sph_hankel > test_sph_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_steed_continued_fraction > test_steed_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_struve > test_struve.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_struve_old > test_struve_old.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_summation > test_summation.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_theta > test_theta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_tr1_cmath > test_tr1_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_tricomi_u > test_tricomi_u.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_trig > test_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_weierstrass_ellint > test_weierstrass_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_wilson > test_wilson.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH ./test_wright_omega > test_wright_omega.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./test_zeta_trig > test_zeta_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./run_coulfg > run_coulfg.txt
	./RUN_COULFG > RUN_COULFG.TXT

check: $(CHECK_DIR) $(CHECKS)
	echo "Beginning executions of checks..." > $(CHECK_DIR)/check_out.txt 2> $(CHECK_DIR)/check_err.txt
	echo "check_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bernoulli" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_bernoulli >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_binomial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_binomial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_u" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_u >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_v" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_v >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_w" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chebyshev_w >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_chi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_clausen_cl" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_clausen_cl >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_digamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_digamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_euler" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_euler >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_eulerian_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_eulerian_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_eulerian_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_eulerian_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_lbinomial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lbinomial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ldouble_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_ldouble_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_legendre_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lfactorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lfactorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_logistic_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_logistic_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lognormal_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lognormal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lfalling_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lfalling_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lrising_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_lrising_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_owens_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_owens_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gamma_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_falling_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_falling_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_rising_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_rising_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gamma_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_stirling_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_stirling_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_stirling_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_stirling_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "pr68397" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/pr68397 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/origin_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


check_tr1: $(CHECK_DIR) $(TR1_CHECKS)
	echo "Beginning executions of tr1 checks..." > $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt
	echo "check_tr1_assoc_laguerre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_assoc_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_beta" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_1" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_comp_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_2" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_comp_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_3" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_comp_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_conf_hyperg" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_i" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_j" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_k" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_cyl_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_neumann" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_1" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_2" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_3" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_expint" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_expint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_hermite" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_hermite >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_hyperg" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_laguerre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_riemann_zeta" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_bessel" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_sph_bessel >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_neumann" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && ${CHECK_DIR}/check_tr1_sph_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


mpfrcalc: mpfr_gexpr.c
	$(GCC) -I. -o mpfrcalc mpfr_gexpr.c -lmpfr -lgmp -lm


test_special_function: wrappers_debug test_special_function.cpp test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -Iwrappers -o test_special_function test_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

diff_special_function: wrappers_debug diff_special_function.cpp test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -Iwrappers -o diff_special_function diff_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase2: wrappers_debug libstdc++_support/testcase2.cpp libstdc++_support/testcase2.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -I. -o libstdc++_support/testcase2 libstdc++_support/testcase2.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -UTR1 -I. -o libstdc++_support/testcase libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase_tr1: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -DTR1 -I. -o libstdc++_support/testcase_tr1 libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

test_Faddeeva: $(FAD_DIR)/Faddeeva.hh $(FAD_DIR)/Faddeeva.cc
	$(CXX) -DTEST_FADDEEVA -o $(FAD_DIR)/test_Faddeeva $(FAD_DIR)/Faddeeva.cc -lquadmath

test_airy_roots: laboratories/airy_functions/test_airy_roots.cpp
	$(CXX17) -I. -o test_airy_roots laboratories/airy_functions/test_airy_roots.cpp -lquadmath

test_anger_weber: laboratories/bessel_functions/test_anger_weber.cpp
	$(CXX17) -I. -o test_anger_weber laboratories/bessel_functions/test_anger_weber.cpp -lquadmath

test_appell_f1: laboratories/appell_functions/test_appell_f1.cpp
	$(CXX17) -I. -o test_appell_f1 laboratories/appell_functions/test_appell_f1.cpp -lquadmath

test_arith_geom_mean: laboratories/norm_functions/test_arith_geom_mean.cpp
	$(CXX17) -I. -Iwrappers -o test_arith_geom_mean laboratories/norm_functions/test_arith_geom_mean.cpp -lquadmath

test_bernoulli: wrappers_debug laboratories/bernoulli_functions/test_bernoulli.cpp
	$(CXX17) -I. -Iwrappers -o test_bernoulli laboratories/bernoulli_functions/test_bernoulli.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

test_bessel: laboratories/bessel_functions/test_bessel.cpp laboratories/bessel_functions/new_bessel.tcc
	$(CXX17) -I. -o test_bessel laboratories/bessel_functions/test_bessel.cpp -lquadmath

test_bessel_asymp: laboratories/bessel_functions/test_bessel_asymp.cpp
	$(CXX17) -I. -o test_bessel_asymp laboratories/bessel_functions/test_bessel_asymp.cpp -lquadmath

test_bessel_iter: laboratories/bessel_functions/test_bessel_iter.cpp
	$(CXX17) -I. -o test_bessel_iter laboratories/bessel_functions/test_bessel_iter.cpp -lquadmath

test_beta: wrappers_debug laboratories/beta_functions/test_beta.cpp
	$(CXX17) -I. -Iwrappers -o test_beta laboratories/beta_functions/test_beta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_beta_inc: laboratories/beta_functions/test_beta_inc.cpp
	$(CXX17) -I. -o test_beta_inc laboratories/beta_functions/test_beta_inc.cpp -lquadmath

test_binet: laboratories/gamma_functions/test_binet.cpp
	$(CXX17) -I. -o test_binet laboratories/gamma_functions/test_binet.cpp -lquadmath

test_binet_float: laboratories/gamma_functions/test_binet_float.cpp
	$(CXX17) -I. -o test_binet_float laboratories/gamma_functions/test_binet_float.cpp -lquadmath

test_bose_einstein: laboratories/zeta_functions/test_bose_einstein.cpp
	$(CXX17) -I. -o test_bose_einstein laboratories/zeta_functions/test_bose_einstein.cpp -lquadmath

test_charlier: laboratories/orthogonal_polynomials/test_charlier.cpp
	$(CXX17) -I. -Iwrappers -o test_charlier laboratories/orthogonal_polynomials/test_charlier.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

test_chebyshev: laboratories/orthogonal_polynomials/test_chebyshev.cpp
	$(CXX17) -I. -o test_chebyshev laboratories/orthogonal_polynomials/test_chebyshev.cpp -lquadmath

test_chebyshev_trig: laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp
	$(CXX17) -I. -o test_chebyshev_trig laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp -lquadmath

test_chebyshev_trig_pi: laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp
	$(CXX17) -I. -o test_chebyshev_trig_pi laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp -lquadmath

test_clausen: wrappers_debug laboratories/zeta_functions/test_clausen.cpp
	$(CXX17) -I. -Iwrappers -o test_clausen laboratories/zeta_functions/test_clausen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_cmath: test_cmath.cpp
	$(CXX) -I. -o test_cmath test_cmath.cpp -lquadmath

test_comp_ellint_1: laboratories/elliptic_integrals/test_comp_ellint_1.cpp
	$(CXX17) -I. -o test_comp_ellint_1 laboratories/elliptic_integrals/test_comp_ellint_1.cpp -lquadmath

test_complex128: test_complex128.cpp
	$(CXX17) -I. -o test_complex128 test_complex128.cpp -lquadmath

test_complex_gamma: laboratories/gamma_functions/test_complex_gamma.cpp
	$(CXX17) -I. -o test_complex_gamma laboratories/gamma_functions/test_complex_gamma.cpp -lquadmath

test_conf_hyperg: laboratories/hypergeometric_functions/test_conf_hyperg.cpp
	$(CXX17) -I. -o test_conf_hyperg laboratories/hypergeometric_functions/test_conf_hyperg.cpp -lquadmath

test_conf_hyperg_limit: laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp
	$(CXX17) -I. -o test_conf_hyperg_limit laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp -lquadmath

test_const: test_const.cpp
	$(CXX17) -I. -I../mpreal -o test_const test_const.cpp -lquadmath -lmpfr -lgmp

test_continued_fraction: test_continued_fraction.cpp
	$(CXX17) -I. -o test_continued_fraction test_continued_fraction.cpp -lquadmath

test_continuous_dual_hahn: laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp
	$(CXX17) -I. -o test_continuous_dual_hahn laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp -lquadmath

test_continuous_hahn: laboratories/orthogonal_polynomials/test_continuous_hahn.cpp
	$(CXX17) -I. -o test_continuous_hahn laboratories/orthogonal_polynomials/test_continuous_hahn.cpp -lquadmath

test_cordic: laboratories/elementary_functions/test_cordic.cpp
	$(CXX17) -I. -o test_cordic laboratories/elementary_functions/test_cordic.cpp -lquadmath

test_coulomb: laboratories/coulomb_functions/test_coulomb.cpp
	$(CXX17) -I. -o test_coulomb laboratories/coulomb_functions/test_coulomb.cpp -lquadmath

test_csint: laboratories/exponential_integrals/test_csint.cpp laboratories/exponential_integrals/csint.tcc
	$(CXX17) -I. -Ilaboratories/exponential_integrals -o test_csint laboratories/exponential_integrals/test_csint.cpp -lquadmath

test_cyl_hankel: wrappers_debug laboratories/bessel_functions/test_cyl_hankel.cpp
	$(CXX17) -I. -Iwrappers -o test_cyl_hankel laboratories/bessel_functions/test_cyl_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_dawson: laboratories/error_functions/test_dawson.cpp
	$(CXX17) -I. -o test_dawson laboratories/error_functions/test_dawson.cpp -lquadmath

test_debye: wrappers_debug laboratories/zeta_functions/test_debye.cpp
	$(CXX17) -I. -Iwrappers -o test_debye laboratories/zeta_functions/test_debye.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_digamma: wrappers_debug laboratories/gamma_functions/test_digamma.cpp
	$(CXX17) -I. -Iwrappers -o test_digamma laboratories/gamma_functions/test_digamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_dilog: laboratories/zeta_functions/test_dilog.cpp
	$(CXX17) -I. -o test_dilog laboratories/zeta_functions/test_dilog.cpp -lquadmath

test_dirichlet_eta: laboratories/zeta_functions/test_dirichlet_eta.cpp
	$(CXX17) -I. -o test_dirichlet_eta laboratories/zeta_functions/test_dirichlet_eta.cpp -lquadmath

test_dual_hahn: laboratories/orthogonal_polynomials/test_dual_hahn.cpp
	$(CXX17) -I. -o test_dual_hahn laboratories/orthogonal_polynomials/test_dual_hahn.cpp -lquadmath

test_erfc: wrappers_debug laboratories/error_functions/test_erfc.cpp
	$(CXX17) -I. -Iwrappers -o test_erfc laboratories/error_functions/test_erfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_experfc: wrappers_debug laboratories/error_functions/test_experfc.cpp
	$(CXX17) -I. -Iwrappers -I../mpreal -o test_experfc laboratories/error_functions/test_experfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost -lmpfr -lgmp

test_expint: wrappers_debug laboratories/exponential_integrals/test_expint.cpp
	$(CXX17) -I. -Iwrappers -o test_expint laboratories/exponential_integrals/test_expint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_factorial: laboratories/gamma_functions/test_factorial.cpp
	$(CXX17) -I. -o test_factorial laboratories/gamma_functions/test_factorial.cpp -lquadmath

test_faddeeva: wrappers_debug laboratories/error_functions/test_faddeeva.cpp
	$(CXX17) -I. -Iwrappers -o test_faddeeva laboratories/error_functions/test_faddeeva.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_faddeeva

test_falling_factorial: wrappers_debug laboratories/gamma_functions/test_falling_factorial.cpp
	$(CXX17) -I. -Iwrappers -o test_falling_factorial laboratories/gamma_functions/test_falling_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_fermi_dirac: wrappers_debug laboratories/zeta_functions/test_fermi_dirac.cpp
	$(CXX17) -I. -Iwrappers -o test_fermi_dirac laboratories/zeta_functions/test_fermi_dirac.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_fibonacci: wrappers_debug laboratories/test_fibonacci.cpp
	$(CXX17) -I. -Iwrappers -o test_fibonacci laboratories/test_fibonacci.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_float128: test_float128.cpp
	$(CXX17) -I. -o test_float128 test_float128.cpp -lquadmath

test_fresnel: wrappers_debug laboratories/error_functions/test_fresnel.cpp laboratories/error_functions/fresnel.tcc
	$(CXX17) -I. -Iwrappers -Ilaboratories/error_functions -o test_fresnel laboratories/error_functions/test_fresnel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_gamma: wrappers_debug laboratories/gamma_functions/test_gamma.cpp
	$(CXX17) -I. -Iwrappers -o test_gamma laboratories/gamma_functions/test_gamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_gamma_ratio: wrappers_debug laboratories/gamma_functions/test_gamma_ratio.cpp
	$(CXX17) -I. -Iwrappers -o test_gamma_ratio laboratories/gamma_functions/test_gamma_ratio.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_gamma_reciprocal: laboratories/gamma_functions/test_gamma_reciprocal.cpp
	$(CXX17) -I. -o test_gamma_reciprocal laboratories/gamma_functions/test_gamma_reciprocal.cpp -lquadmath

test_gegenbauer: laboratories/orthogonal_polynomials/test_gegenbauer.cpp
	$(CXX17) -I. -o test_gegenbauer laboratories/orthogonal_polynomials/test_gegenbauer.cpp -lquadmath

test_gudermannian: laboratories/elementary_functions/test_gudermannian.cpp
	$(CXX17) -I. -o test_gudermannian laboratories/elementary_functions/test_gudermannian.cpp -lquadmath

test_hahn: laboratories/orthogonal_polynomials/test_hahn.cpp
	$(CXX17) -I. -o test_hahn laboratories/orthogonal_polynomials/test_hahn.cpp -lquadmath

test_hankel: laboratories/bessel_functions/test_hankel.cpp
	$(CXX17) -I. -o test_hankel laboratories/bessel_functions/test_hankel.cpp -lquadmath

test_hankel_real_arg: wrappers_debug laboratories/bessel_functions/test_hankel_real_arg.cpp
	$(CXX17) -I. -Iwrappers -o test_hankel_real_arg laboratories/bessel_functions/test_hankel_real_arg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_hermite: laboratories/orthogonal_polynomials/test_hermite.cpp laboratories/orthogonal_polynomials/new_hermite.tcc
	$(CXX17) -I. -Ilaboratories/orthogonal_polynomials -o test_hermite laboratories/orthogonal_polynomials/test_hermite.cpp -lquadmath

test_heuman_lambda: wrappers_debug laboratories/elliptic_integrals/test_heuman_lambda.cpp
	$(CXX17) -I. -Iwrappers -o test_heuman_lambda laboratories/elliptic_integrals/test_heuman_lambda.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_hurwitz_zeta: laboratories/zeta_functions/test_hurwitz_zeta.cpp
	$(CXX17) -I. -o test_hurwitz_zeta laboratories/zeta_functions/test_hurwitz_zeta.cpp -lquadmath

test_hurwitz_zeta_new: laboratories/zeta_functions/test_hurwitz_zeta_new.cpp
	$(CXX17) -I. -o test_hurwitz_zeta_new laboratories/zeta_functions/test_hurwitz_zeta_new.cpp -lquadmath

test_hydrogen: wrappers_debug laboratories/coulomb_functions/test_hydrogen.cpp
	$(CXX17) -I. -Iwrappers -o test_hydrogen laboratories/coulomb_functions/test_hydrogen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_hyperg: wrappers_debug laboratories/hypergeometric_functions/test_hyperg.cpp
	$(CXX17) -I. -Iwrappers -o test_hyperg laboratories/hypergeometric_functions/test_hyperg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_hypot: laboratories/norm_functions/test_hypot.cpp
	$(CXX17) -I. -o test_hypot laboratories/norm_functions/test_hypot.cpp -lquadmath

test_inv_erf: test_inv_erf.cpp
	$(CXX17) -I. -DSTANDALONE -o test_inv_erf test_inv_erf.cpp -lquadmath

test_inv_gamma: laboratories/gamma_functions/test_inv_gamma.cpp test_inv_erf.cpp
	$(CXX17) -I. -o test_inv_gamma laboratories/gamma_functions/test_inv_gamma.cpp -lquadmath

test_inv_ibeta: laboratories/beta_functions/test_inv_ibeta.cpp
	$(CXX17) -I. -o test_inv_ibeta laboratories/beta_functions/test_inv_ibeta.cpp -lquadmath

test_jacobi: laboratories/orthogonal_polynomials/test_jacobi.cpp
	$(CXX17) -I. -o test_jacobi laboratories/orthogonal_polynomials/test_jacobi.cpp -lquadmath

test_jacobi_ellint: laboratories/theta_functions/test_jacobi_ellint.cpp
	$(CXX17) -I. -Iwrappers -o test_jacobi_ellint laboratories/theta_functions/test_jacobi_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

test_jacobi_inv: laboratories/theta_functions/test_jacobi_inv.cpp
	$(CXX17) -I. -o test_jacobi_inv laboratories/theta_functions/test_jacobi_inv.cpp -lquadmath

test_jacobi_theta: wrappers_debug laboratories/theta_functions/test_jacobi_theta.cpp
	$(CXX17) -I. -Iwrappers -o test_jacobi_theta laboratories/theta_functions/test_jacobi_theta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_jacobi_zeta: wrappers_debug laboratories/elliptic_integrals/test_jacobi_zeta.cpp
	$(CXX17) -I. -Iwrappers -o test_jacobi_zeta laboratories/elliptic_integrals/test_jacobi_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_kelvin: laboratories/bessel_functions/test_kelvin.cpp
	$(CXX17) -I. -o test_kelvin laboratories/bessel_functions/test_kelvin.cpp -lquadmath

test_krawtchouk: laboratories/orthogonal_polynomials/test_krawtchouk.cpp
	$(CXX17) -I. -Iwrappers -o test_krawtchouk laboratories/orthogonal_polynomials/test_krawtchouk.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

test_laguerre: laboratories/orthogonal_polynomials/test_laguerre.cpp
	$(CXX17) -I. -o test_laguerre laboratories/orthogonal_polynomials/test_laguerre.cpp -lquadmath

test_lambert_w: laboratories/elementary_functions/test_lambert_w.cpp
	$(CXX17) -I. -o test_lambert_w laboratories/elementary_functions/test_lambert_w.cpp -lquadmath

test_large_order_bessel: laboratories/bessel_functions/test_large_order_bessel.cpp
	$(CXX17) -I. -o test_large_order_bessel laboratories/bessel_functions/test_large_order_bessel.cpp -lquadmath

test_legendre: laboratories/orthogonal_polynomials/test_legendre.cpp
	$(CXX17) -I. -o test_legendre laboratories/orthogonal_polynomials/test_legendre.cpp -lquadmath

test_legendre_ellint: laboratories/elliptic_integrals/test_legendre_ellint.cpp
	$(CXX17) -I. -Iwrappers -o test_legendre_ellint laboratories/elliptic_integrals/test_legendre_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_lentz_continued_fraction: test_lentz_continued_fraction.cpp
	$(CXX17) -I. -o test_lentz_continued_fraction test_lentz_continued_fraction.cpp -lquadmath

test_lerch: laboratories/zeta_functions/test_lerch.cpp
	$(CXX17) -I. -o test_lerch laboratories/zeta_functions/test_lerch.cpp lerchphi/Source/lerchphi.cpp -lquadmath

test_limits: test_limits.cpp
	$(CXX17) -I. -o test_limits test_limits.cpp -lquadmath

test_little_airy: wrappers_debug laboratories/airy_functions/test_little_airy.cpp
	$(CXX17) -I. -Iwrappers -o test_little_airy laboratories/airy_functions/test_little_airy.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_lobatto: wrappers_debug laboratories/orthogonal_polynomials/test_lobatto.cpp
	$(CXX17) -I. -Iwrappers -o test_lobatto laboratories/orthogonal_polynomials/test_lobatto.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_lommel: wrappers_debug laboratories/bessel_functions/test_lommel.cpp
	$(CXX17) -I. -Iwrappers -o test_lommel laboratories/bessel_functions/test_lommel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_logsumexp: laboratories/norm_functions/test_logsumexp.cpp
	$(CXX17) -I. -o test_logsumexp laboratories/norm_functions/test_logsumexp.cpp -lquadmath

test_marcum_q: laboratories/distributions/test_marcum_q.cpp
	$(CXX17) -I. -o test_marcum_q laboratories/distributions/test_marcum_q.cpp -lquadmath

test_math_h: test_math_h.cpp
	$(CXX) -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o test_math_h test_math_h.cpp -lquadmath

test_maxint: laboratories/floating_point_tools/test_maxint.cpp
	$(CXX17) -I. -I../mpreal -o test_maxint laboratories/floating_point_tools/test_maxint.cpp -lquadmath -lmpfr

test_meixner: laboratories/orthogonal_polynomials/test_meixner.cpp
	$(CXX17) -I. -Iwrappers -o test_meixner laboratories/orthogonal_polynomials/test_meixner.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

test_meixner_pollaczek: laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp
	$(CXX17) -I. -Iwrappers -o test_meixner_pollaczek laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp -lquadmath

test_mittag_leffler: laboratories/mittag_leffler_functions/test_mittag_leffler.cpp
	$(CXX17) -I. -o test_mittag_leffler laboratories/mittag_leffler_functions/test_mittag_leffler.cpp -lquadmath

test_mod2pi: test_mod2pi.cpp
	$(CXX17) -I. -I../mpreal -o test_mod2pi test_mod2pi.cpp -lquadmath -lmpfr -lgmp

test_mpreal: test_mpreal.cpp
	$(CXX17) -I. -I../mpreal -o test_mpreal test_mpreal.cpp -lquadmath -lmpfr -lgmp

test_notsospecfun: test_notsospecfun.cpp
	$(CXX17) -I. -o test_notsospecfun test_notsospecfun.cpp -lquadmath

test_nric_bessel: laboratories/bessel_functions/test_nric_bessel.cpp laboratories/bessel_functions/nric_bessel.tcc
	$(CXX) -I. -Ilaboratories/bessel_functions -o test_nric_bessel laboratories/bessel_functions/test_nric_bessel.cpp -lquadmath

test_numeric_limits: test_numeric_limits.cpp
	$(CXX17) -I. -I../mpreal -o test_numeric_limits test_numeric_limits.cpp -lquadmath -lmpfr -lgmp

test_owens_t: wrappers_debug laboratories/error_functions/test_owens_t.cpp
	$(CXX17) -I. -Iwrappers -o test_owens_t laboratories/error_functions/test_owens_t.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_parab_cyl: laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp
	$(CXX17) -I. -o test_parab_cyl laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp -lquadmath

test_polygamma: laboratories/gamma_functions/test_polygamma.cpp
	$(CXX17) -I. -Iwrappers -o test_polygamma laboratories/gamma_functions/test_polygamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_polylog: wrappers_debug laboratories/zeta_functions/test_polylog.cpp
	$(CXX17) -I. -Iwrappers -o test_polylog laboratories/zeta_functions/test_polylog.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_cephes

test_power_mean: laboratories/norm_functions/test_power_mean.cpp
	$(CXX17) -I. -o test_power_mean laboratories/norm_functions/test_power_mean.cpp -lquadmath

test_power_norm: laboratories/norm_functions/test_power_norm.cpp
	$(CXX17) -I. -o test_power_norm laboratories/norm_functions/test_power_norm.cpp -lquadmath

test_racah: laboratories/orthogonal_polynomials/test_racah.cpp
	$(CXX17) -I. -o test_racah laboratories/orthogonal_polynomials/test_racah.cpp -lquadmath

test_rational: test_rational.cpp
	$(CXX17) -I. -o test_rational test_rational.cpp -lquadmath

test_recursion: test_recursion.cpp
	$(CXX17) -I. -o test_recursion test_recursion.cpp -lquadmath

test_reperiodized_hyper: laboratories/elementary_functions/test_reperiodized_hyper.cpp
	$(CXX17) -I. -o test_reperiodized_hyper laboratories/elementary_functions/test_reperiodized_hyper.cpp -lquadmath

test_reperiodized_trig: wrappers_debug laboratories/elementary_functions/test_reperiodized_trig.cpp
	$(CXX17) -I. -Iwrappers -o test_reperiodized_trig laboratories/elementary_functions/test_reperiodized_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_riemann_zeta: laboratories/zeta_functions/test_riemann_zeta.cpp
	$(CXX17) -I. -Iwrappers -o test_riemann_zeta laboratories/zeta_functions/test_riemann_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

test_rising_factorial: wrappers_debug laboratories/gamma_functions/test_rising_factorial.cpp
	$(CXX17) -I. -Iwrappers -o test_rising_factorial laboratories/gamma_functions/test_rising_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_root_finding: test_root_finding.cpp
	$(CXX17) -I. -o test_root_finding test_root_finding.cpp -lquadmath

test_sincos: laboratories/elementary_functions/test_sincos.cpp
	$(CXX17) -I. -o test_sincos laboratories/elementary_functions/test_sincos.cpp -lquadmath

test_sinus_cardinal: wrappers_debug laboratories/elementary_functions/test_sinus_cardinal.cpp
	$(CXX17) -I. -Iwrappers -o test_sinus_cardinal laboratories/elementary_functions/test_sinus_cardinal.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

test_sph_bessel: laboratories/bessel_functions/test_sph_bessel.cpp
	$(CXX17) -I. -o test_sph_bessel laboratories/bessel_functions/test_sph_bessel.cpp -lquadmath

test_sph_hankel: wrappers_debug laboratories/bessel_functions/test_sph_hankel.cpp
	$(CXX17) -I. -Iwrappers -o test_sph_hankel laboratories/bessel_functions/test_sph_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

test_steed_continued_fraction: test_steed_continued_fraction.cpp
	$(CXX17) -I. -o test_steed_continued_fraction test_steed_continued_fraction.cpp -lquadmath

test_struve: laboratories/bessel_functions/test_struve.cpp
	$(CXX17) -I. -Iwrappers -o test_struve laboratories/bessel_functions/test_struve.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

test_struve_old: laboratories/bessel_functions/test_struve_old.cpp
	$(CXX17) -I. -Ilaboratories/hypergeometric_functions -o test_struve_old laboratories/bessel_functions/test_struve_old.cpp -lquadmath

test_summation: test_summation.cpp
	$(CXX17) -I. -o test_summation test_summation.cpp -lquadmath

test_theta: laboratories/theta_functions/test_theta.cpp
	$(CXX17) -I. -o test_theta laboratories/theta_functions/test_theta.cpp -lquadmath

test_tr1_cmath: test_tr1_cmath.cpp
	$(CXX) -I. -o test_tr1_cmath test_tr1_cmath.cpp -lquadmath

test_tricomi_u: laboratories/hypergeometric_functions/test_tricomi_u.cpp
	$(CXX17) -I. -o test_tricomi_u laboratories/hypergeometric_functions/test_tricomi_u.cpp -lquadmath

test_trig: laboratories/elementary_functions/test_trig.cpp
	$(CXX17) -I. -o test_trig laboratories/elementary_functions/test_trig.cpp -lquadmath

test_weierstrass_ellint: laboratories/theta_functions/test_weierstrass_ellint.cpp
	$(CXX17) -I. -o test_weierstrass_ellint laboratories/theta_functions/test_weierstrass_ellint.cpp -lquadmath

test_wilson: laboratories/orthogonal_polynomials/test_wilson.cpp
	$(CXX17) -I. -o test_wilson laboratories/orthogonal_polynomials/test_wilson.cpp -lquadmath

test_wright_omega: laboratories/elementary_functions/test_wright_omega.cpp
	$(CXX17) -I. -o test_wright_omega laboratories/elementary_functions/test_wright_omega.cpp -lquadmath

test_zeta_trig: laboratories/elementary_functions/test_zeta_trig.cpp
	$(CXX17) -I. -o test_zeta_trig laboratories/elementary_functions/test_zeta_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl


run_coulfg: laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp
	$(CXX17) -I. -o run_coulfg laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp -lquadmath

RUN_COULFG: laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR
	gfortran -o RUN_COULFG laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR


airy_toy: laboratories/airy_functions/airy_toy.cpp
	$(CXX17) -I. -o airy_toy laboratories/airy_functions/airy_toy.cpp -lquadmath

hankel_toy: laboratories/bessel_functions/hankel_toy.cpp
	$(CXX17) -I. -o hankel_toy laboratories/bessel_functions/hankel_toy.cpp -lquadmath

hankel_toy128: laboratories/bessel_functions/hankel_toy128.cpp
	$(CXX17) -I. -o hankel_toy128 laboratories/bessel_functions/hankel_toy128.cpp -lquadmath

hankel_toy_new: laboratories/bessel_functions/hankel_toy_new.cpp
	$(CXX17) -I. -o hankel_toy_new laboratories/bessel_functions/hankel_toy_new.cpp -lquadmath


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

${CHECK_DIR}/check_binomial: ${CHECK_DIR}/check_binomial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_binomial ${CHECK_DIR}/check_binomial.cc -lquadmath

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

${CHECK_DIR}/check_clausen_cl: ${CHECK_DIR}/check_clausen_cl.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_clausen_cl ${CHECK_DIR}/check_clausen_cl.cc -lquadmath

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

${CHECK_DIR}/check_digamma: ${CHECK_DIR}/check_digamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_digamma ${CHECK_DIR}/check_digamma.cc -lquadmath

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

${CHECK_DIR}/check_euler: ${CHECK_DIR}/check_euler.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_euler ${CHECK_DIR}/check_euler.cc -lquadmath

${CHECK_DIR}/check_eulerian_1: ${CHECK_DIR}/check_eulerian_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_eulerian_1 ${CHECK_DIR}/check_eulerian_1.cc -lquadmath

${CHECK_DIR}/check_eulerian_2: ${CHECK_DIR}/check_eulerian_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_eulerian_2 ${CHECK_DIR}/check_eulerian_2.cc -lquadmath

${CHECK_DIR}/check_expint: ${CHECK_DIR}/check_expint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint ${CHECK_DIR}/check_expint.cc -lquadmath

${CHECK_DIR}/check_expint_en: ${CHECK_DIR}/check_expint_en.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint_en ${CHECK_DIR}/check_expint_en.cc -lquadmath

${CHECK_DIR}/check_factorial: ${CHECK_DIR}/check_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_factorial ${CHECK_DIR}/check_factorial.cc -lquadmath

${CHECK_DIR}/check_falling_factorial: ${CHECK_DIR}/check_falling_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_falling_factorial ${CHECK_DIR}/check_falling_factorial.cc -lquadmath

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

${CHECK_DIR}/check_lbinomial: ${CHECK_DIR}/check_lbinomial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lbinomial ${CHECK_DIR}/check_lbinomial.cc -lquadmath

${CHECK_DIR}/check_ldouble_factorial: ${CHECK_DIR}/check_ldouble_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ldouble_factorial ${CHECK_DIR}/check_ldouble_factorial.cc -lquadmath

${CHECK_DIR}/check_legendre: ${CHECK_DIR}/check_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre ${CHECK_DIR}/check_legendre.cc -lquadmath

${CHECK_DIR}/check_legendre_q: ${CHECK_DIR}/check_legendre_q.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre_q ${CHECK_DIR}/check_legendre_q.cc -lquadmath

${CHECK_DIR}/check_lfactorial: ${CHECK_DIR}/check_lfactorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfactorial ${CHECK_DIR}/check_lfactorial.cc -lquadmath

${CHECK_DIR}/check_lfalling_factorial: ${CHECK_DIR}/check_lfalling_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfalling_factorial ${CHECK_DIR}/check_lfalling_factorial.cc -lquadmath

${CHECK_DIR}/check_lgamma: ${CHECK_DIR}/check_lgamma.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lgamma ${CHECK_DIR}/check_lgamma.cc -lquadmath

${CHECK_DIR}/check_logistic_p: ${CHECK_DIR}/check_logistic_p.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_p ${CHECK_DIR}/check_logistic_p.cc -lquadmath

${CHECK_DIR}/check_logistic_pdf: ${CHECK_DIR}/check_logistic_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_pdf ${CHECK_DIR}/check_logistic_pdf.cc -lquadmath

${CHECK_DIR}/check_lognormal_p: ${CHECK_DIR}/check_lognormal_p.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_p ${CHECK_DIR}/check_lognormal_p.cc -lquadmath

${CHECK_DIR}/check_lognormal_pdf: ${CHECK_DIR}/check_lognormal_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_pdf ${CHECK_DIR}/check_lognormal_pdf.cc -lquadmath

${CHECK_DIR}/check_lrising_factorial: ${CHECK_DIR}/check_lrising_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lrising_factorial ${CHECK_DIR}/check_lrising_factorial.cc -lquadmath

${CHECK_DIR}/check_normal_p: ${CHECK_DIR}/check_normal_p.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_p ${CHECK_DIR}/check_normal_p.cc -lquadmath

${CHECK_DIR}/check_normal_pdf: ${CHECK_DIR}/check_normal_pdf.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_pdf ${CHECK_DIR}/check_normal_pdf.cc -lquadmath

${CHECK_DIR}/check_owens_t: ${CHECK_DIR}/check_owens_t.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_owens_t ${CHECK_DIR}/check_owens_t.cc -lquadmath

${CHECK_DIR}/check_gamma_p: ${CHECK_DIR}/check_gamma_p.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_p ${CHECK_DIR}/check_gamma_p.cc -lquadmath

${CHECK_DIR}/check_gamma_q: ${CHECK_DIR}/check_gamma_q.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_q ${CHECK_DIR}/check_gamma_q.cc -lquadmath

${CHECK_DIR}/check_radpoly: ${CHECK_DIR}/check_radpoly.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_radpoly ${CHECK_DIR}/check_radpoly.cc -lquadmath

${CHECK_DIR}/check_riemann_zeta: ${CHECK_DIR}/check_riemann_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_riemann_zeta ${CHECK_DIR}/check_riemann_zeta.cc -lquadmath

${CHECK_DIR}/check_rising_factorial: ${CHECK_DIR}/check_rising_factorial.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_rising_factorial ${CHECK_DIR}/check_rising_factorial.cc -lquadmath

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

${CHECK_DIR}/check_stirling_1: ${CHECK_DIR}/check_stirling_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_stirling_1 ${CHECK_DIR}/check_stirling_1.cc -lquadmath

${CHECK_DIR}/check_stirling_2: ${CHECK_DIR}/check_stirling_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_stirling_2 ${CHECK_DIR}/check_stirling_2.cc -lquadmath

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

${CHECK_DIR}/pr68397: ${CHECK_DIR}/pr68397.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr68397 ${CHECK_DIR}/pr68397.cc -lquadmath

${CHECK_DIR}/origin_cyl_bessel_j: ${CHECK_DIR}/origin_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_bessel_j ${CHECK_DIR}/origin_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/origin_cyl_neumann: ${CHECK_DIR}/origin_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_neumann ${CHECK_DIR}/origin_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_tr1_assoc_laguerre: ${CHECK_DIR}/check_tr1_assoc_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_assoc_laguerre ${CHECK_DIR}/check_tr1_assoc_laguerre.cc -lquadmath

${CHECK_DIR}/check_tr1_assoc_legendre: ${CHECK_DIR}/check_tr1_assoc_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_assoc_legendre ${CHECK_DIR}/check_tr1_assoc_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_beta: ${CHECK_DIR}/check_tr1_beta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_beta ${CHECK_DIR}/check_tr1_beta.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_1: ${CHECK_DIR}/check_tr1_comp_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_1 ${CHECK_DIR}/check_tr1_comp_ellint_1.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_2: ${CHECK_DIR}/check_tr1_comp_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_2 ${CHECK_DIR}/check_tr1_comp_ellint_2.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_3: ${CHECK_DIR}/check_tr1_comp_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_3 ${CHECK_DIR}/check_tr1_comp_ellint_3.cc -lquadmath

${CHECK_DIR}/check_tr1_conf_hyperg: ${CHECK_DIR}/check_tr1_conf_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_conf_hyperg ${CHECK_DIR}/check_tr1_conf_hyperg.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_i: ${CHECK_DIR}/check_tr1_cyl_bessel_i.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_i ${CHECK_DIR}/check_tr1_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_j: ${CHECK_DIR}/check_tr1_cyl_bessel_j.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_j ${CHECK_DIR}/check_tr1_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_k: ${CHECK_DIR}/check_tr1_cyl_bessel_k.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_k ${CHECK_DIR}/check_tr1_cyl_bessel_k.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_neumann: ${CHECK_DIR}/check_tr1_cyl_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_neumann ${CHECK_DIR}/check_tr1_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_1: ${CHECK_DIR}/check_tr1_ellint_1.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_1 ${CHECK_DIR}/check_tr1_ellint_1.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_2: ${CHECK_DIR}/check_tr1_ellint_2.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_2 ${CHECK_DIR}/check_tr1_ellint_2.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_3: ${CHECK_DIR}/check_tr1_ellint_3.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_3 ${CHECK_DIR}/check_tr1_ellint_3.cc -lquadmath

${CHECK_DIR}/check_tr1_expint: ${CHECK_DIR}/check_tr1_expint.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_expint ${CHECK_DIR}/check_tr1_expint.cc -lquadmath

${CHECK_DIR}/check_tr1_hermite: ${CHECK_DIR}/check_tr1_hermite.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_hermite ${CHECK_DIR}/check_tr1_hermite.cc -lquadmath

${CHECK_DIR}/check_tr1_hyperg: ${CHECK_DIR}/check_tr1_hyperg.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_hyperg ${CHECK_DIR}/check_tr1_hyperg.cc -lquadmath

${CHECK_DIR}/check_tr1_laguerre: ${CHECK_DIR}/check_tr1_laguerre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_laguerre ${CHECK_DIR}/check_tr1_laguerre.cc -lquadmath

${CHECK_DIR}/check_tr1_legendre: ${CHECK_DIR}/check_tr1_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_legendre ${CHECK_DIR}/check_tr1_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_riemann_zeta: ${CHECK_DIR}/check_tr1_riemann_zeta.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_riemann_zeta ${CHECK_DIR}/check_tr1_riemann_zeta.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_bessel: ${CHECK_DIR}/check_tr1_sph_bessel.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_bessel ${CHECK_DIR}/check_tr1_sph_bessel.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_legendre: ${CHECK_DIR}/check_tr1_sph_legendre.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_legendre ${CHECK_DIR}/check_tr1_sph_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_neumann: ${CHECK_DIR}/check_tr1_sph_neumann.cc
	$(CXX) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_neumann ${CHECK_DIR}/check_tr1_sph_neumann.cc -lquadmath


$(CHECK_DIR): $(CHECK_DIR)
	if test ! -d $(CHECK_DIR); then \
	  mkdir $(CHECK_DIR); \
	fi


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
