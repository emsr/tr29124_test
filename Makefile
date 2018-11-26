
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
CXX14 = $(CXX_INST_DIR)/bin/g++ -std=gnu++14 -g -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wall -Wextra -Wno-psabi
CXX17 = $(CXX_INST_DIR)/bin/g++ -std=gnu++17 -fconcepts -g -Wall -Wextra -Wno-psabi
CXX2A = $(CXX_INST_DIR)/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = libstdc++_support

INC_DIR = include/bits
INCLUDES = -Iinclude -Ipolynomial/include

FAD_DIR = 3rdparty/Faddeeva

MATH_DIR = ./libstdc++_support

CHECK_DIR = $(MATH_DIR)/check

OBJ_DIR = $(MATH_DIR)/obj

LIB_DIR = $(MATH_DIR)

TEST_BIN_DIR = laboratories/bin
TEST_OUT_DIR = laboratories/output
TEST_SF_OUT_DIR = laboratories/test_output
DIFF_SF_OUT_DIR = laboratories/diff_output


#       test_special_function \
#       testcase2 \
#

BINS = \
       libstdc++_support/testcase \
       libstdc++_support/testcase_tr1 \
       $(TEST_BIN_DIR)/mpfrcalc \
       $(TEST_BIN_DIR)/diff_special_function \
       $(TEST_BIN_DIR)/test_special_function \
       $(TEST_BIN_DIR)/airy_toy \
       $(TEST_BIN_DIR)/hankel_toy \
       $(TEST_BIN_DIR)/hankel_toy128 \
       $(TEST_BIN_DIR)/hankel_toy_new \
       $(TEST_BIN_DIR)/test_airy_roots \
       $(TEST_BIN_DIR)/test_anger_weber \
       $(TEST_BIN_DIR)/test_appell_f1 \
       $(TEST_BIN_DIR)/test_arith_geom_mean \
       $(TEST_BIN_DIR)/test_bernoulli \
       $(TEST_BIN_DIR)/test_bessel \
       $(TEST_BIN_DIR)/test_bessel_asymp \
       $(TEST_BIN_DIR)/test_bessel_iter \
       $(TEST_BIN_DIR)/test_beta \
       $(TEST_BIN_DIR)/test_beta_inc \
       $(TEST_BIN_DIR)/test_binet \
       $(TEST_BIN_DIR)/test_binet_float \
       $(TEST_BIN_DIR)/test_bose_einstein \
       $(TEST_BIN_DIR)/test_charlier \
       $(TEST_BIN_DIR)/test_chebyshev \
       $(TEST_BIN_DIR)/test_chebyshev_trig \
       $(TEST_BIN_DIR)/test_chebyshev_trig_pi \
       $(TEST_BIN_DIR)/test_clausen \
       $(TEST_BIN_DIR)/test_cmath \
       $(TEST_BIN_DIR)/test_comp_ellint_1 \
       $(TEST_BIN_DIR)/test_complex128 \
       $(TEST_BIN_DIR)/test_complex_gamma \
       $(TEST_BIN_DIR)/test_conf_hyperg \
       $(TEST_BIN_DIR)/test_conf_hyperg_limit \
       $(TEST_BIN_DIR)/test_const \
       $(TEST_BIN_DIR)/test_continued_fraction \
       $(TEST_BIN_DIR)/test_continuous_dual_hahn \
       $(TEST_BIN_DIR)/test_continuous_hahn \
       $(TEST_BIN_DIR)/test_cordic \
       $(TEST_BIN_DIR)/test_coulomb \
       $(TEST_BIN_DIR)/test_csint \
       $(TEST_BIN_DIR)/test_cyl_hankel \
       $(TEST_BIN_DIR)/test_dawson \
       $(TEST_BIN_DIR)/test_debye \
       $(TEST_BIN_DIR)/test_digamma \
       $(TEST_BIN_DIR)/test_dilog \
       $(TEST_BIN_DIR)/test_dirichlet_eta \
       $(TEST_BIN_DIR)/test_dual_hahn \
       $(TEST_BIN_DIR)/test_erfc \
       $(TEST_BIN_DIR)/test_experfc \
       $(TEST_BIN_DIR)/test_expint \
       $(TEST_BIN_DIR)/test_factorial \
       $(TEST_BIN_DIR)/test_faddeeva \
       $(TEST_BIN_DIR)/test_falling_factorial \
       $(TEST_BIN_DIR)/test_fermi_dirac \
       $(TEST_BIN_DIR)/test_fibonacci \
       $(TEST_BIN_DIR)/test_float128 \
       $(TEST_BIN_DIR)/test_fresnel \
       $(TEST_BIN_DIR)/test_gamma \
       $(TEST_BIN_DIR)/test_gamma_ratio \
       $(TEST_BIN_DIR)/test_gamma_reciprocal \
       $(TEST_BIN_DIR)/test_gegenbauer \
       $(TEST_BIN_DIR)/test_gudermannian \
       $(TEST_BIN_DIR)/test_hahn \
       $(TEST_BIN_DIR)/test_hankel \
       $(TEST_BIN_DIR)/test_hankel_real_arg \
       $(TEST_BIN_DIR)/test_hermite \
       $(TEST_BIN_DIR)/test_heuman_lambda \
       $(TEST_BIN_DIR)/test_hurwitz_zeta \
       $(TEST_BIN_DIR)/test_hurwitz_zeta_new \
       $(TEST_BIN_DIR)/test_hydrogen \
       $(TEST_BIN_DIR)/test_hyperg \
       $(TEST_BIN_DIR)/test_hypot \
       $(TEST_BIN_DIR)/test_inv_erf \
       $(TEST_BIN_DIR)/test_inv_gamma \
       $(TEST_BIN_DIR)/test_inv_ibeta \
       $(TEST_BIN_DIR)/test_inv_lgamma \
       $(TEST_BIN_DIR)/test_jacobi \
       $(TEST_BIN_DIR)/test_jacobi_ellint \
       $(TEST_BIN_DIR)/test_jacobi_inv \
       $(TEST_BIN_DIR)/test_jacobi_theta \
       $(TEST_BIN_DIR)/test_jacobi_zeta \
       $(TEST_BIN_DIR)/test_kelvin \
       $(TEST_BIN_DIR)/test_krawtchouk \
       $(TEST_BIN_DIR)/test_laguerre \
       $(TEST_BIN_DIR)/test_lambert_w \
       $(TEST_BIN_DIR)/test_large_order_bessel \
       $(TEST_BIN_DIR)/test_legendre \
       $(TEST_BIN_DIR)/test_legendre_ellint \
       $(TEST_BIN_DIR)/test_lentz_continued_fraction \
       $(TEST_BIN_DIR)/test_lerch \
       $(TEST_BIN_DIR)/test_limits \
       $(TEST_BIN_DIR)/test_little_airy \
       $(TEST_BIN_DIR)/test_lobatto \
       $(TEST_BIN_DIR)/test_logsumexp \
       $(TEST_BIN_DIR)/test_lommel \
       $(TEST_BIN_DIR)/test_marcum_q \
       $(TEST_BIN_DIR)/test_math_h \
       $(TEST_BIN_DIR)/test_maxint \
       $(TEST_BIN_DIR)/test_meixner \
       $(TEST_BIN_DIR)/test_meixner_pollaczek \
       $(TEST_BIN_DIR)/test_mittag_leffler \
       $(TEST_BIN_DIR)/test_mod2pi \
       $(TEST_BIN_DIR)/test_mpreal \
       $(TEST_BIN_DIR)/test_notsospecfun \
       $(TEST_BIN_DIR)/test_nric_bessel \
       $(TEST_BIN_DIR)/test_numeric_limits \
       $(TEST_BIN_DIR)/test_owens_t \
       $(TEST_BIN_DIR)/test_parab_cyl \
       $(TEST_BIN_DIR)/test_polygamma \
       $(TEST_BIN_DIR)/test_polylog \
       $(TEST_BIN_DIR)/test_power_mean \
       $(TEST_BIN_DIR)/test_power_norm \
       $(TEST_BIN_DIR)/test_pow_limits \
       $(TEST_BIN_DIR)/test_racah \
       $(TEST_BIN_DIR)/test_rational \
       $(TEST_BIN_DIR)/test_recursion \
       $(TEST_BIN_DIR)/test_reperiodized_hyper \
       $(TEST_BIN_DIR)/test_reperiodized_trig \
       $(TEST_BIN_DIR)/test_riemann_zeta \
       $(TEST_BIN_DIR)/test_rising_factorial \
       $(TEST_BIN_DIR)/test_root_search \
       $(TEST_BIN_DIR)/test_sincos \
       $(TEST_BIN_DIR)/test_sinus_cardinal \
       $(TEST_BIN_DIR)/test_sph_bessel \
       $(TEST_BIN_DIR)/test_sph_hankel \
       $(TEST_BIN_DIR)/test_steed_continued_fraction \
       $(TEST_BIN_DIR)/test_struve \
       $(TEST_BIN_DIR)/test_struve_old \
       $(TEST_BIN_DIR)/test_summation \
       $(TEST_BIN_DIR)/test_theta \
       $(TEST_BIN_DIR)/test_tr1_cmath \
       $(TEST_BIN_DIR)/test_tricomi_u \
       $(TEST_BIN_DIR)/test_trig \
       $(TEST_BIN_DIR)/test_weierstrass_ellint \
       $(TEST_BIN_DIR)/test_wilson \
       $(TEST_BIN_DIR)/test_wright_omega \
       $(TEST_BIN_DIR)/test_zeta_trig \
       $(TEST_BIN_DIR)/run_coulfg \
       $(TEST_BIN_DIR)/RUN_COULFG


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
	 ${CHECK_DIR}/check_clausen_cl \
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
	 ${CHECK_DIR}/check_debye \
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
	 ${CHECK_DIR}/check_falling_factorial \
	 ${CHECK_DIR}/check_fresnel_c \
	 ${CHECK_DIR}/check_fresnel_s \
	 ${CHECK_DIR}/check_gamma_p \
	 ${CHECK_DIR}/check_gamma_q \
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
	 ${CHECK_DIR}/check_lfalling_factorial \
	 ${CHECK_DIR}/check_lgamma \
	 ${CHECK_DIR}/check_logistic_p \
	 ${CHECK_DIR}/check_logistic_pdf \
	 ${CHECK_DIR}/check_lrising_factorial \
	 ${CHECK_DIR}/check_lognormal_p \
	 ${CHECK_DIR}/check_lognormal_pdf \
	 ${CHECK_DIR}/check_normal_p \
	 ${CHECK_DIR}/check_normal_pdf \
	 ${CHECK_DIR}/check_owens_t \
	 ${CHECK_DIR}/check_polygamma \
	 ${CHECK_DIR}/check_radpoly \
	 ${CHECK_DIR}/check_riemann_zeta \
	 ${CHECK_DIR}/check_rising_factorial \
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


all: $(TEST_BIN_DIR) $(BINS)


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


docs: include/bits/*
	rm -rf docs/html/*
	rm -rf docs/latex/*
	doxygen
	cd docs/latex && make

testcases2: libstdc++_support/testcase2
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./libstdc++_support/testcase2

testcases: libstdc++_support/testcase
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./libstdc++_support/testcase

testcases_tr1: libstdc++_support/testcase_tr1
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./libstdc++_support/testcase_tr1

diffs: $(DIFF_SF_OUT_DIR) laboratories/bin/diff_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./laboratories/bin/diff_special_function > $(DIFF_SF_OUT_DIR)/diff_special_function.txt

tests: $(TEST_SF_OUT_DIR) laboratories/bin/test_special_function
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH ./laboratories/bin/test_special_function > $(TEST_SF_OUT_DIR)/test_special_function.txt

# This will always build the executables!
test: $(BINS) $(TEST_OUT_DIR)
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/airy_toy > $(TEST_OUT_DIR)/airy_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy > $(TEST_OUT_DIR)/hankel_toy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy128 > $(TEST_OUT_DIR)/hankel_toy128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy_new > $(TEST_OUT_DIR)/hankel_toy_new.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_airy_roots > $(TEST_OUT_DIR)/test_airy_roots.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_anger_weber > $(TEST_OUT_DIR)/test_anger_weber.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_appell_f1 > $(TEST_OUT_DIR)/test_appell_f1.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bernoulli > $(TEST_OUT_DIR)/test_bernoulli.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel > $(TEST_OUT_DIR)/test_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel_asymp > $(TEST_OUT_DIR)/test_bessel_asymp.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel_iter > $(TEST_OUT_DIR)/test_bessel_iter.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_beta > $(TEST_OUT_DIR)/test_beta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_beta_inc > $(TEST_OUT_DIR)/test_beta_inc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_binet > $(TEST_OUT_DIR)/test_binet.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_binet_float > $(TEST_OUT_DIR)/test_binet_float.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bose_einstein > $(TEST_OUT_DIR)/test_bose_einstein.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_charlier > $(TEST_OUT_DIR)/test_charlier.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev > $(TEST_OUT_DIR)/test_chebyshev.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev_trig > $(TEST_OUT_DIR)/test_chebyshev_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev_trig_pi > $(TEST_OUT_DIR)/test_chebyshev_trig_pi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_clausen > $(TEST_OUT_DIR)/test_clausen.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cmath > $(TEST_OUT_DIR)/test_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_comp_ellint_1 > $(TEST_OUT_DIR)/test_comp_ellint_1.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_complex128 > $(TEST_OUT_DIR)/test_complex128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_complex_gamma > $(TEST_OUT_DIR)/test_complex_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_conf_hyperg > $(TEST_OUT_DIR)/test_conf_hyperg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_conf_hyperg_limit > $(TEST_OUT_DIR)/test_conf_hyperg_limit.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_const > $(TEST_OUT_DIR)/test_const.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continued_fraction > $(TEST_OUT_DIR)/test_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continuous_dual_hahn > $(TEST_OUT_DIR)/test_continuous_dual_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continuous_hahn > $(TEST_OUT_DIR)/test_continuous_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cordic > $(TEST_OUT_DIR)/test_cordic.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_coulomb > $(TEST_OUT_DIR)/test_coulomb.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_csint > $(TEST_OUT_DIR)/test_csint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cyl_hankel > $(TEST_OUT_DIR)/test_cyl_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dawson > $(TEST_OUT_DIR)/test_dawson.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_debye > $(TEST_OUT_DIR)/test_debye.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_digamma > $(TEST_OUT_DIR)/test_digamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dilog > $(TEST_OUT_DIR)/test_dilog.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dirichlet_eta > $(TEST_OUT_DIR)/test_dirichlet_eta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dual_hahn > $(TEST_OUT_DIR)/test_dual_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_erfc > $(TEST_OUT_DIR)/test_erfc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_experfc > $(TEST_OUT_DIR)/test_experfc.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_expint > $(TEST_OUT_DIR)/test_expint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_factorial > $(TEST_OUT_DIR)/test_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_faddeeva > $(TEST_OUT_DIR)/test_faddeeva.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_falling_factorial > $(TEST_OUT_DIR)/test_falling_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fermi_dirac > $(TEST_OUT_DIR)/test_fermi_dirac.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fibonacci > $(TEST_OUT_DIR)/test_fibonacci.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_float128 > $(TEST_OUT_DIR)/test_float128.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fresnel > $(TEST_OUT_DIR)/test_fresnel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma > $(TEST_OUT_DIR)/test_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma_ratio > $(TEST_OUT_DIR)/test_gamma_ratio.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma_reciprocal > $(TEST_OUT_DIR)/test_gamma_reciprocal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gegenbauer > $(TEST_OUT_DIR)/test_gegenbauer.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gudermannian > $(TEST_OUT_DIR)/test_gudermannian.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hahn > $(TEST_OUT_DIR)/test_hahn.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hankel > $(TEST_OUT_DIR)/test_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hankel_real_arg > $(TEST_OUT_DIR)/test_hankel_real_arg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hermite > $(TEST_OUT_DIR)/test_hermite.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_heuman_lambda > $(TEST_OUT_DIR)/test_heuman_lambda.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hurwitz_zeta > $(TEST_OUT_DIR)/test_hurwitz_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hurwitz_zeta_new > $(TEST_OUT_DIR)/test_hurwitz_zeta_new.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hydrogen > $(TEST_OUT_DIR)/test_hydrogen.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hyperg > $(TEST_OUT_DIR)/test_hyperg.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hypot > $(TEST_OUT_DIR)/test_hypot.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_erf > $(TEST_OUT_DIR)/test_inv_erf.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_gamma > $(TEST_OUT_DIR)/test_inv_gamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_ibeta > $(TEST_OUT_DIR)/test_inv_ibeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_lgamma > $(TEST_OUT_DIR)/test_inv_lgamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi > $(TEST_OUT_DIR)/test_jacobi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_ellint > $(TEST_OUT_DIR)/test_jacobi_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_inv > $(TEST_OUT_DIR)/test_jacobi_inv.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_theta > $(TEST_OUT_DIR)/test_jacobi_theta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_zeta > $(TEST_OUT_DIR)/test_jacobi_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_kelvin > $(TEST_OUT_DIR)/test_kelvin.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_krawtchouk > $(TEST_OUT_DIR)/test_krawtchouk.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_laguerre > $(TEST_OUT_DIR)/test_laguerre.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lambert_w > $(TEST_OUT_DIR)/test_lambert_w.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_legendre > $(TEST_OUT_DIR)/test_legendre.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_legendre_ellint > $(TEST_OUT_DIR)/test_legendre_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lentz_continued_fraction > $(TEST_OUT_DIR)/test_lentz_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lerch > $(TEST_OUT_DIR)/test_lerch.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_limits > $(TEST_OUT_DIR)/test_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_little_airy > $(TEST_OUT_DIR)/test_little_airy.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lobatto > $(TEST_OUT_DIR)/test_lobatto.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_logsumexp > $(TEST_OUT_DIR)/test_logsumexp.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lommel > $(TEST_OUT_DIR)/test_lommel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_marcum_q > $(TEST_OUT_DIR)/test_marcum_q.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_math_h > $(TEST_OUT_DIR)/test_math_h.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_maxint > $(TEST_OUT_DIR)/test_maxint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_meixner > $(TEST_OUT_DIR)/test_meixner.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_meixner_pollaczek > $(TEST_OUT_DIR)/test_meixner_pollaczek.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mittag_leffler > $(TEST_OUT_DIR)/test_mittag_leffler.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mod2pi > $(TEST_OUT_DIR)/test_mod2pi.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mpreal > $(TEST_OUT_DIR)/test_mpreal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_notsospecfun > $(TEST_OUT_DIR)/test_notsospecfun.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_nric_bessel > $(TEST_OUT_DIR)/test_nric_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_numeric_limits > $(TEST_OUT_DIR)/test_numeric_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_owens_t > $(TEST_OUT_DIR)/test_owens_t.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_parab_cyl > $(TEST_OUT_DIR)/test_parab_cyl.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_polygamma > $(TEST_OUT_DIR)/test_polygamma.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_polylog > $(TEST_OUT_DIR)/test_polylog.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_power_mean > $(TEST_OUT_DIR)/test_power_mean.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_power_norm > $(TEST_OUT_DIR)/test_power_norm.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_pow_limits > $(TEST_OUT_DIR)/test_pow_limits.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_racah > $(TEST_OUT_DIR)/test_racah.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_rational > $(TEST_OUT_DIR)/test_rational.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_recursion > $(TEST_OUT_DIR)/test_recursion.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_reperiodized_hyper > $(TEST_OUT_DIR)/test_reperiodized_hyper.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_reperiodized_trig > $(TEST_OUT_DIR)/test_reperiodized_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_riemann_zeta > $(TEST_OUT_DIR)/test_riemann_zeta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_rising_factorial > $(TEST_OUT_DIR)/test_rising_factorial.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_root_search > $(TEST_OUT_DIR)/test_root_search.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sincos > $(TEST_OUT_DIR)/test_sincos.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sinus_cardinal > $(TEST_OUT_DIR)/test_sinus_cardinal.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sph_bessel > $(TEST_OUT_DIR)/test_sph_bessel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sph_hankel > $(TEST_OUT_DIR)/test_sph_hankel.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_steed_continued_fraction > $(TEST_OUT_DIR)/test_steed_continued_fraction.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_struve > $(TEST_OUT_DIR)/test_struve.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_struve_old > $(TEST_OUT_DIR)/test_struve_old.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_summation > $(TEST_OUT_DIR)/test_summation.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_theta > $(TEST_OUT_DIR)/test_theta.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_tr1_cmath > $(TEST_OUT_DIR)/test_tr1_cmath.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_tricomi_u > $(TEST_OUT_DIR)/test_tricomi_u.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_trig > $(TEST_OUT_DIR)/test_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_weierstrass_ellint > $(TEST_OUT_DIR)/test_weierstrass_ellint.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_wilson > $(TEST_OUT_DIR)/test_wilson.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_wright_omega > $(TEST_OUT_DIR)/test_wright_omega.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_zeta_trig > $(TEST_OUT_DIR)/test_zeta_trig.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/run_coulfg > $(TEST_OUT_DIR)/run_coulfg.txt
	./RUN_COULFG > $(TEST_OUT_DIR)/RUN_COULFG.TXT

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
	echo "check_debye" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_debye >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_falling_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_falling_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_c" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_c >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_s" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_fresnel_s >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gamma_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_gamma_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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
	echo "check_normal_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_normal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_owens_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_owens_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_polygamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_polygamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_radpoly" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_radpoly >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_riemann_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_rising_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && $(CHECK_DIR)/check_rising_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
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


libstdc++_support/testcase2: wrappers_debug libstdc++_support/testcase2.cpp libstdc++_support/testcase2.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) $(INCLUDES) -o libstdc++_support/testcase2 libstdc++_support/testcase2.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -UTR1 $(INCLUDES) -o libstdc++_support/testcase libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase_tr1: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) -DTR1 $(INCLUDES) -o libstdc++_support/testcase_tr1 libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva


$(TEST_BIN_DIR)/mpfrcalc: multiprecision/mpfr_gexpr.c
	$(GCC) -Iinclude -o $(TEST_BIN_DIR)/mpfrcalc multiprecision/mpfr_gexpr.c -lmpfr -lgmp -lm

$(TEST_BIN_DIR)/test_special_function: wrappers_debug laboratories/test_special_function.cpp laboratories/test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_special_function laboratories/test_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

$(TEST_BIN_DIR)/diff_special_function: wrappers_debug laboratories/diff_special_function.cpp laboratories/test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/diff_special_function laboratories/diff_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

$(TEST_BIN_DIR)/test_Faddeeva: $(FAD_DIR)/Faddeeva.hh $(FAD_DIR)/Faddeeva.cc
	$(CXX14) -DTEST_FADDEEVA -o $(TEST_BIN_DIR)/$(FAD_DIR)/test_Faddeeva $(FAD_DIR)/Faddeeva.cc -lquadmath

$(TEST_BIN_DIR)/test_airy_roots: laboratories/airy_functions/test_airy_roots.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_airy_roots laboratories/airy_functions/test_airy_roots.cpp -lquadmath

$(TEST_BIN_DIR)/test_anger_weber: laboratories/bessel_functions/test_anger_weber.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_anger_weber laboratories/bessel_functions/test_anger_weber.cpp -lquadmath

$(TEST_BIN_DIR)/test_appell_f1: laboratories/appell_functions/test_appell_f1.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_appell_f1 laboratories/appell_functions/test_appell_f1.cpp -lquadmath

$(TEST_BIN_DIR)/test_arith_geom_mean: laboratories/norm_functions/test_arith_geom_mean.cpp
	$(CXX17) -Iinclude -Iwrappers -o $(TEST_BIN_DIR)/test_arith_geom_mean laboratories/norm_functions/test_arith_geom_mean.cpp -lquadmath

$(TEST_BIN_DIR)/test_bernoulli: wrappers_debug laboratories/bernoulli_functions/test_bernoulli.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_bernoulli laboratories/bernoulli_functions/test_bernoulli.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

$(TEST_BIN_DIR)/test_bessel: laboratories/bessel_functions/test_bessel.cpp laboratories/bessel_functions/new_bessel.tcc
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel laboratories/bessel_functions/test_bessel.cpp -lquadmath

$(TEST_BIN_DIR)/test_bessel_asymp: laboratories/bessel_functions/test_bessel_asymp.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel_asymp laboratories/bessel_functions/test_bessel_asymp.cpp -lquadmath

$(TEST_BIN_DIR)/test_bessel_iter: laboratories/bessel_functions/test_bessel_iter.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel_iter laboratories/bessel_functions/test_bessel_iter.cpp -lquadmath

$(TEST_BIN_DIR)/test_beta: wrappers_debug laboratories/beta_functions/test_beta.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_beta laboratories/beta_functions/test_beta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_beta_inc: laboratories/beta_functions/test_beta_inc.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_beta_inc laboratories/beta_functions/test_beta_inc.cpp -lquadmath

$(TEST_BIN_DIR)/test_binet: laboratories/gamma_functions/test_binet.cpp
	$(CXX17) $(INCLUDES) -I. -Irational/include -o $(TEST_BIN_DIR)/test_binet laboratories/gamma_functions/test_binet.cpp -lquadmath

$(TEST_BIN_DIR)/test_binet_float: laboratories/gamma_functions/test_binet_float.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_binet_float laboratories/gamma_functions/test_binet_float.cpp -lquadmath

$(TEST_BIN_DIR)/test_bose_einstein: laboratories/zeta_functions/test_bose_einstein.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bose_einstein laboratories/zeta_functions/test_bose_einstein.cpp -lquadmath

$(TEST_BIN_DIR)/test_charlier: laboratories/orthogonal_polynomials/test_charlier.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_charlier laboratories/orthogonal_polynomials/test_charlier.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

$(TEST_BIN_DIR)/test_chebyshev: laboratories/orthogonal_polynomials/test_chebyshev.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev laboratories/orthogonal_polynomials/test_chebyshev.cpp -lquadmath

$(TEST_BIN_DIR)/test_chebyshev_trig: laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev_trig laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp -lquadmath

$(TEST_BIN_DIR)/test_chebyshev_trig_pi: laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev_trig_pi laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp -lquadmath

$(TEST_BIN_DIR)/test_clausen: wrappers_debug laboratories/zeta_functions/test_clausen.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_clausen laboratories/zeta_functions/test_clausen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_cmath: test_std_maths/test_cmath.cpp
	$(CXX14) $(INCLUDES) -o $(TEST_BIN_DIR)/test_cmath test_std_maths/test_cmath.cpp -lquadmath

$(TEST_BIN_DIR)/test_comp_ellint_1: laboratories/elliptic_integrals/test_comp_ellint_1.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_comp_ellint_1 laboratories/elliptic_integrals/test_comp_ellint_1.cpp -lquadmath

$(TEST_BIN_DIR)/test_complex128: laboratories/complex_tools/test_complex128.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_complex128 laboratories/complex_tools/test_complex128.cpp -lquadmath

$(TEST_BIN_DIR)/test_complex_gamma: laboratories/gamma_functions/test_complex_gamma.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_complex_gamma laboratories/gamma_functions/test_complex_gamma.cpp -lquadmath

$(TEST_BIN_DIR)/test_conf_hyperg: laboratories/hypergeometric_functions/test_conf_hyperg.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_conf_hyperg laboratories/hypergeometric_functions/test_conf_hyperg.cpp -lquadmath

$(TEST_BIN_DIR)/test_conf_hyperg_limit: laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_conf_hyperg_limit laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp -lquadmath

$(TEST_BIN_DIR)/test_const: laboratories/constants/test_const.cpp
	$(CXX17) $(INCLUDES) -I../mpreal -o $(TEST_BIN_DIR)/test_const laboratories/constants/test_const.cpp -lquadmath -lmpfr -lgmp

$(TEST_BIN_DIR)/test_continued_fraction: continued_fractions/test_continued_fraction.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continued_fraction continued_fractions/test_continued_fraction.cpp -lquadmath

$(TEST_BIN_DIR)/test_continuous_dual_hahn: laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continuous_dual_hahn laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp -lquadmath

$(TEST_BIN_DIR)/test_continuous_hahn: laboratories/orthogonal_polynomials/test_continuous_hahn.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continuous_hahn laboratories/orthogonal_polynomials/test_continuous_hahn.cpp -lquadmath

$(TEST_BIN_DIR)/test_cordic: laboratories/elementary_functions/test_cordic.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_cordic laboratories/elementary_functions/test_cordic.cpp -lquadmath

$(TEST_BIN_DIR)/test_coulomb: laboratories/coulomb_functions/test_coulomb.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_coulomb laboratories/coulomb_functions/test_coulomb.cpp -lquadmath

$(TEST_BIN_DIR)/test_csint: laboratories/exponential_integrals/test_csint.cpp laboratories/exponential_integrals/csint.tcc
	$(CXX17) $(INCLUDES) -Ilaboratories/exponential_integrals -o $(TEST_BIN_DIR)/test_csint laboratories/exponential_integrals/test_csint.cpp -lquadmath

$(TEST_BIN_DIR)/test_cyl_hankel: wrappers_debug laboratories/bessel_functions/test_cyl_hankel.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_cyl_hankel laboratories/bessel_functions/test_cyl_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_dawson: laboratories/error_functions/test_dawson.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dawson laboratories/error_functions/test_dawson.cpp -lquadmath

$(TEST_BIN_DIR)/test_debye: wrappers_debug laboratories/zeta_functions/test_debye.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_debye laboratories/zeta_functions/test_debye.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_digamma: wrappers_debug laboratories/gamma_functions/test_digamma.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_digamma laboratories/gamma_functions/test_digamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_dilog: laboratories/zeta_functions/test_dilog.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dilog laboratories/zeta_functions/test_dilog.cpp -lquadmath

$(TEST_BIN_DIR)/test_dirichlet_eta: laboratories/zeta_functions/test_dirichlet_eta.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dirichlet_eta laboratories/zeta_functions/test_dirichlet_eta.cpp -lquadmath

$(TEST_BIN_DIR)/test_dual_hahn: laboratories/orthogonal_polynomials/test_dual_hahn.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dual_hahn laboratories/orthogonal_polynomials/test_dual_hahn.cpp -lquadmath

$(TEST_BIN_DIR)/test_erfc: wrappers_debug laboratories/error_functions/test_erfc.cpp
	$(CXX17) $(INCLUDES) -Icontinued_fractions/include -Iwrappers -o $(TEST_BIN_DIR)/test_erfc laboratories/error_functions/test_erfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_experfc: wrappers_debug laboratories/error_functions/test_experfc.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -I../mpreal -o $(TEST_BIN_DIR)/test_experfc laboratories/error_functions/test_experfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost -lmpfr -lgmp

$(TEST_BIN_DIR)/test_expint: wrappers_debug laboratories/exponential_integrals/test_expint.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_expint laboratories/exponential_integrals/test_expint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_factorial: laboratories/gamma_functions/test_factorial.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_factorial laboratories/gamma_functions/test_factorial.cpp -lquadmath

$(TEST_BIN_DIR)/test_faddeeva: wrappers_debug laboratories/error_functions/test_faddeeva.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_faddeeva laboratories/error_functions/test_faddeeva.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_faddeeva

$(TEST_BIN_DIR)/test_falling_factorial: wrappers_debug laboratories/gamma_functions/test_falling_factorial.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_falling_factorial laboratories/gamma_functions/test_falling_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_fermi_dirac: wrappers_debug laboratories/zeta_functions/test_fermi_dirac.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_fermi_dirac laboratories/zeta_functions/test_fermi_dirac.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_fibonacci: wrappers_debug laboratories/test_fibonacci.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_fibonacci laboratories/test_fibonacci.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_float128: laboratories/floating_point_tools/test_float128.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_float128 laboratories/floating_point_tools/test_float128.cpp -lquadmath

$(TEST_BIN_DIR)/test_fresnel: wrappers_debug laboratories/error_functions/test_fresnel.cpp laboratories/error_functions/fresnel.tcc
	$(CXX17) $(INCLUDES) -Iwrappers -Ilaboratories/error_functions -o $(TEST_BIN_DIR)/test_fresnel laboratories/error_functions/test_fresnel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_gamma: wrappers_debug laboratories/gamma_functions/test_gamma.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_gamma laboratories/gamma_functions/test_gamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_gamma_ratio: wrappers_debug laboratories/gamma_functions/test_gamma_ratio.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_gamma_ratio laboratories/gamma_functions/test_gamma_ratio.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_gamma_reciprocal: laboratories/gamma_functions/test_gamma_reciprocal.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_gamma_reciprocal laboratories/gamma_functions/test_gamma_reciprocal.cpp -lquadmath

$(TEST_BIN_DIR)/test_gegenbauer: laboratories/orthogonal_polynomials/test_gegenbauer.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_gegenbauer laboratories/orthogonal_polynomials/test_gegenbauer.cpp -lquadmath

$(TEST_BIN_DIR)/test_gudermannian: laboratories/elementary_functions/test_gudermannian.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_gudermannian laboratories/elementary_functions/test_gudermannian.cpp -lquadmath

$(TEST_BIN_DIR)/test_hahn: laboratories/orthogonal_polynomials/test_hahn.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hahn laboratories/orthogonal_polynomials/test_hahn.cpp -lquadmath

$(TEST_BIN_DIR)/test_hankel: laboratories/bessel_functions/test_hankel.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hankel laboratories/bessel_functions/test_hankel.cpp -lquadmath

$(TEST_BIN_DIR)/test_hankel_real_arg: wrappers_debug laboratories/bessel_functions/test_hankel_real_arg.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hankel_real_arg laboratories/bessel_functions/test_hankel_real_arg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_hermite: laboratories/orthogonal_polynomials/test_hermite.cpp laboratories/orthogonal_polynomials/new_hermite.tcc
	$(CXX17) $(INCLUDES) -Iquadrature/include -Icontinued_fractions/include -Ilaboratories/orthogonal_polynomials -o $(TEST_BIN_DIR)/test_hermite laboratories/orthogonal_polynomials/test_hermite.cpp -lquadmath

$(TEST_BIN_DIR)/test_heuman_lambda: wrappers_debug laboratories/elliptic_integrals/test_heuman_lambda.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_heuman_lambda laboratories/elliptic_integrals/test_heuman_lambda.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_hurwitz_zeta: laboratories/zeta_functions/test_hurwitz_zeta.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hurwitz_zeta laboratories/zeta_functions/test_hurwitz_zeta.cpp -lquadmath

$(TEST_BIN_DIR)/test_hurwitz_zeta_new: laboratories/zeta_functions/test_hurwitz_zeta_new.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hurwitz_zeta_new laboratories/zeta_functions/test_hurwitz_zeta_new.cpp -lquadmath

$(TEST_BIN_DIR)/test_hydrogen: wrappers_debug laboratories/coulomb_functions/test_hydrogen.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hydrogen laboratories/coulomb_functions/test_hydrogen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_hyperg: wrappers_debug laboratories/hypergeometric_functions/test_hyperg.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hyperg laboratories/hypergeometric_functions/test_hyperg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_hypot: laboratories/norm_functions/test_hypot.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hypot laboratories/norm_functions/test_hypot.cpp -lquadmath

$(TEST_BIN_DIR)/test_inv_erf: laboratories/error_functions/test_inv_erf.cpp
	$(CXX17) $(INCLUDES) -DSTANDALONE -o $(TEST_BIN_DIR)/test_inv_erf laboratories/error_functions/test_inv_erf.cpp -lquadmath

$(TEST_BIN_DIR)/test_inv_gamma: laboratories/gamma_functions/test_inv_gamma.cpp laboratories/error_functions/test_inv_erf.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_inv_gamma laboratories/gamma_functions/test_inv_gamma.cpp -lquadmath

$(TEST_BIN_DIR)/test_inv_ibeta: laboratories/beta_functions/test_inv_ibeta.cpp
	$(CXX17) $(INCLUDES) -Iroot_search/include -o $(TEST_BIN_DIR)/test_inv_ibeta laboratories/beta_functions/test_inv_ibeta.cpp -lquadmath

$(TEST_BIN_DIR)/test_inv_lgamma: laboratories/gamma_functions/test_inv_lgamma.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_inv_lgamma laboratories/gamma_functions/test_inv_lgamma.cpp -lquadmath

$(TEST_BIN_DIR)/test_jacobi: laboratories/orthogonal_polynomials/test_jacobi.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_jacobi laboratories/orthogonal_polynomials/test_jacobi.cpp -lquadmath

$(TEST_BIN_DIR)/test_jacobi_ellint: laboratories/theta_functions/test_jacobi_ellint.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_ellint laboratories/theta_functions/test_jacobi_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

$(TEST_BIN_DIR)/test_jacobi_inv: laboratories/theta_functions/test_jacobi_inv.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_jacobi_inv laboratories/theta_functions/test_jacobi_inv.cpp -lquadmath

$(TEST_BIN_DIR)/test_jacobi_theta: wrappers_debug laboratories/theta_functions/test_jacobi_theta.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_theta laboratories/theta_functions/test_jacobi_theta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_jacobi_zeta: wrappers_debug laboratories/elliptic_integrals/test_jacobi_zeta.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_zeta laboratories/elliptic_integrals/test_jacobi_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_kelvin: laboratories/bessel_functions/test_kelvin.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_kelvin laboratories/bessel_functions/test_kelvin.cpp -lquadmath

$(TEST_BIN_DIR)/test_krawtchouk: laboratories/orthogonal_polynomials/test_krawtchouk.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_krawtchouk laboratories/orthogonal_polynomials/test_krawtchouk.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

$(TEST_BIN_DIR)/test_laguerre: laboratories/orthogonal_polynomials/test_laguerre.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_laguerre laboratories/orthogonal_polynomials/test_laguerre.cpp -lquadmath

$(TEST_BIN_DIR)/test_lambert_w: laboratories/elementary_functions/test_lambert_w.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_lambert_w laboratories/elementary_functions/test_lambert_w.cpp -lquadmath

$(TEST_BIN_DIR)/test_large_order_bessel: laboratories/bessel_functions/test_large_order_bessel.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_large_order_bessel laboratories/bessel_functions/test_large_order_bessel.cpp -lquadmath

$(TEST_BIN_DIR)/test_legendre: laboratories/orthogonal_polynomials/test_legendre.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_legendre laboratories/orthogonal_polynomials/test_legendre.cpp -lquadmath

$(TEST_BIN_DIR)/test_legendre_ellint: laboratories/elliptic_integrals/test_legendre_ellint.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_legendre_ellint laboratories/elliptic_integrals/test_legendre_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_lentz_continued_fraction: continued_fractions/test_lentz_continued_fraction.cpp
	$(CXX17) $(INCLUDES) -Icontinued_fractions/include -o $(TEST_BIN_DIR)/test_lentz_continued_fraction continued_fractions/test_lentz_continued_fraction.cpp -lquadmath

$(TEST_BIN_DIR)/test_lerch: laboratories/zeta_functions/test_lerch.cpp
	$(CXX17) $(INCLUDES) -I. -o $(TEST_BIN_DIR)/test_lerch laboratories/zeta_functions/test_lerch.cpp 3rdparty/lerchphi/Source/lerchphi.cpp -lquadmath

$(TEST_BIN_DIR)/test_limits: laboratories/floating_point_tools/test_limits.cpp
	$(CXX17) -Iinclude -o $(TEST_BIN_DIR)/test_limits laboratories/floating_point_tools/test_limits.cpp -lquadmath

$(TEST_BIN_DIR)/test_little_airy: wrappers_debug laboratories/airy_functions/test_little_airy.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_little_airy laboratories/airy_functions/test_little_airy.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_lobatto: wrappers_debug laboratories/orthogonal_polynomials/test_lobatto.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_lobatto laboratories/orthogonal_polynomials/test_lobatto.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_lommel: wrappers_debug laboratories/bessel_functions/test_lommel.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_lommel laboratories/bessel_functions/test_lommel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_logsumexp: laboratories/norm_functions/test_logsumexp.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_logsumexp laboratories/norm_functions/test_logsumexp.cpp -lquadmath

$(TEST_BIN_DIR)/test_marcum_q: laboratories/distributions/test_marcum_q.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_marcum_q laboratories/distributions/test_marcum_q.cpp -lquadmath

$(TEST_BIN_DIR)/test_math_h: test_std_maths/test_math_h.cpp
	$(CXX14) $(INCLUDES) -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o $(TEST_BIN_DIR)/test_math_h test_std_maths/test_math_h.cpp -lquadmath

$(TEST_BIN_DIR)/test_maxint: laboratories/floating_point_tools/test_maxint.cpp
	$(CXX17) $(INCLUDES) -I../mpreal -o $(TEST_BIN_DIR)/test_maxint laboratories/floating_point_tools/test_maxint.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/test_meixner: laboratories/orthogonal_polynomials/test_meixner.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_meixner laboratories/orthogonal_polynomials/test_meixner.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

$(TEST_BIN_DIR)/test_meixner_pollaczek: laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_meixner_pollaczek laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp -lquadmath

$(TEST_BIN_DIR)/test_mittag_leffler: laboratories/mittag_leffler_functions/test_mittag_leffler.cpp
	$(CXX17) $(INCLUDES) -Iquadrature/include -o $(TEST_BIN_DIR)/test_mittag_leffler laboratories/mittag_leffler_functions/test_mittag_leffler.cpp -lquadmath

$(TEST_BIN_DIR)/test_mod2pi: laboratories/floating_point_tools/test_mod2pi.cpp
	$(CXX17) $(INCLUDES) -I../mpreal -o $(TEST_BIN_DIR)/test_mod2pi laboratories/floating_point_tools/test_mod2pi.cpp -lquadmath -lmpfr -lgmp

$(TEST_BIN_DIR)/test_mpreal: multiprecision/test_mpreal.cpp
	$(CXX17) $(INCLUDES) -I../mpreal -o $(TEST_BIN_DIR)/test_mpreal multiprecision/test_mpreal.cpp -lquadmath -lmpfr -lgmp

$(TEST_BIN_DIR)/test_notsospecfun: laboratories/elementary_functions/test_notsospecfun.cpp
	$(CXX17) -Iinclude -o $(TEST_BIN_DIR)/test_notsospecfun laboratories/elementary_functions/test_notsospecfun.cpp -lquadmath

$(TEST_BIN_DIR)/test_nric_bessel: laboratories/bessel_functions/test_nric_bessel.cpp laboratories/bessel_functions/nric_bessel.tcc
	$(CXX14) $(INCLUDES) -Ilaboratories/bessel_functions -o $(TEST_BIN_DIR)/test_nric_bessel laboratories/bessel_functions/test_nric_bessel.cpp -lquadmath

$(TEST_BIN_DIR)/test_numeric_limits: laboratories/floating_point_tools/test_numeric_limits.cpp
	$(CXX17) $(INCLUDES) -I../mpreal -o $(TEST_BIN_DIR)/test_numeric_limits laboratories/floating_point_tools/test_numeric_limits.cpp -lquadmath -lmpfr -lgmp

$(TEST_BIN_DIR)/test_owens_t: wrappers_debug laboratories/error_functions/test_owens_t.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_owens_t laboratories/error_functions/test_owens_t.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_parab_cyl: laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_parab_cyl laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp -lquadmath

$(TEST_BIN_DIR)/test_polygamma: laboratories/gamma_functions/test_polygamma.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -Icontinued_fractions/include -I. -o $(TEST_BIN_DIR)/test_polygamma laboratories/gamma_functions/test_polygamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_polylog: wrappers_debug laboratories/zeta_functions/test_polylog.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_polylog laboratories/zeta_functions/test_polylog.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_cephes

$(TEST_BIN_DIR)/test_power_mean: laboratories/norm_functions/test_power_mean.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_power_mean laboratories/norm_functions/test_power_mean.cpp -lquadmath

$(TEST_BIN_DIR)/test_power_norm: laboratories/norm_functions/test_power_norm.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_power_norm laboratories/norm_functions/test_power_norm.cpp -lquadmath

$(TEST_BIN_DIR)/test_pow_limits: test_std_maths/test_pow_limits.cpp
	$(CXX17) -DBIT -o $(TEST_BIN_DIR)/test_pow_limits test_std_maths/test_pow_limits.cpp -lquadmath

$(TEST_BIN_DIR)/test_racah: laboratories/orthogonal_polynomials/test_racah.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_racah laboratories/orthogonal_polynomials/test_racah.cpp -lquadmath

$(TEST_BIN_DIR)/test_rational: rational/test_rational.cpp
	$(CXX17) -Irational/include $(INCLUDES) -o $(TEST_BIN_DIR)/test_rational rational/test_rational.cpp -lquadmath

$(TEST_BIN_DIR)/test_recursion: recursion/test_recursion.cpp
	$(CXX17) -Iinclude -o $(TEST_BIN_DIR)/test_recursion recursion/test_recursion.cpp -lquadmath

$(TEST_BIN_DIR)/test_reperiodized_hyper: laboratories/elementary_functions/test_reperiodized_hyper.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_reperiodized_hyper laboratories/elementary_functions/test_reperiodized_hyper.cpp -lquadmath

$(TEST_BIN_DIR)/test_reperiodized_trig: wrappers_debug laboratories/elementary_functions/test_reperiodized_trig.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_reperiodized_trig laboratories/elementary_functions/test_reperiodized_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_riemann_zeta: laboratories/zeta_functions/test_riemann_zeta.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_riemann_zeta laboratories/zeta_functions/test_riemann_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

$(TEST_BIN_DIR)/test_rising_factorial: wrappers_debug laboratories/gamma_functions/test_rising_factorial.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_rising_factorial laboratories/gamma_functions/test_rising_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_root_search: root_search/test_root_search.cpp root_search/include/ext/*
	$(CXX17) $(INCLUDES) -Iroot_search/include -o $(TEST_BIN_DIR)/test_root_search root_search/test_root_search.cpp -lquadmath

$(TEST_BIN_DIR)/test_sincos: laboratories/elementary_functions/test_sincos.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_sincos laboratories/elementary_functions/test_sincos.cpp -lquadmath

$(TEST_BIN_DIR)/test_sinus_cardinal: wrappers_debug laboratories/elementary_functions/test_sinus_cardinal.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_sinus_cardinal laboratories/elementary_functions/test_sinus_cardinal.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

$(TEST_BIN_DIR)/test_sph_bessel: laboratories/bessel_functions/test_sph_bessel.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_sph_bessel laboratories/bessel_functions/test_sph_bessel.cpp -lquadmath

$(TEST_BIN_DIR)/test_sph_hankel: wrappers_debug laboratories/bessel_functions/test_sph_hankel.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_sph_hankel laboratories/bessel_functions/test_sph_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

$(TEST_BIN_DIR)/test_steed_continued_fraction: continued_fractions/test_steed_continued_fraction.cpp
	$(CXX17) $(INCLUDES) -Icontinued_fractions/include -o $(TEST_BIN_DIR)/test_steed_continued_fraction continued_fractions/test_steed_continued_fraction.cpp -lquadmath

$(TEST_BIN_DIR)/test_struve: laboratories/bessel_functions/test_struve.cpp
	$(CXX17) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_struve laboratories/bessel_functions/test_struve.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

$(TEST_BIN_DIR)/test_struve_old: laboratories/bessel_functions/test_struve_old.cpp
	$(CXX17) $(INCLUDES) -Ilaboratories/hypergeometric_functions -o $(TEST_BIN_DIR)/test_struve_old laboratories/bessel_functions/test_struve_old.cpp -lquadmath

$(TEST_BIN_DIR)/test_summation: summation/test_summation.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_summation summation/test_summation.cpp -lquadmath

$(TEST_BIN_DIR)/test_theta: laboratories/theta_functions/test_theta.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_theta laboratories/theta_functions/test_theta.cpp -lquadmath

$(TEST_BIN_DIR)/test_tr1_cmath: test_std_maths/test_tr1_cmath.cpp
	$(CXX14) $(INCLUDES) -o $(TEST_BIN_DIR)/test_tr1_cmath test_std_maths/test_tr1_cmath.cpp -lquadmath

$(TEST_BIN_DIR)/test_tricomi_u: laboratories/hypergeometric_functions/test_tricomi_u.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_tricomi_u laboratories/hypergeometric_functions/test_tricomi_u.cpp -lquadmath

$(TEST_BIN_DIR)/test_trig: laboratories/elementary_functions/test_trig.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_trig laboratories/elementary_functions/test_trig.cpp -lquadmath

$(TEST_BIN_DIR)/test_weierstrass_ellint: laboratories/theta_functions/test_weierstrass_ellint.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_weierstrass_ellint laboratories/theta_functions/test_weierstrass_ellint.cpp -lquadmath

$(TEST_BIN_DIR)/test_wilson: laboratories/orthogonal_polynomials/test_wilson.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_wilson laboratories/orthogonal_polynomials/test_wilson.cpp -lquadmath

$(TEST_BIN_DIR)/test_wright_omega: laboratories/elementary_functions/test_wright_omega.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_wright_omega laboratories/elementary_functions/test_wright_omega.cpp -lquadmath

$(TEST_BIN_DIR)/test_zeta_trig: laboratories/elementary_functions/test_zeta_trig.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/test_zeta_trig laboratories/elementary_functions/test_zeta_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl


$(TEST_BIN_DIR)/run_coulfg: laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/run_coulfg laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp -lquadmath

$(TEST_BIN_DIR)/RUN_COULFG: laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR
	gfortran -o $(TEST_BIN_DIR)/RUN_COULFG laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR


$(TEST_BIN_DIR)/airy_toy: laboratories/airy_functions/airy_toy.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/airy_toy laboratories/airy_functions/airy_toy.cpp -lquadmath

$(TEST_BIN_DIR)/hankel_toy: laboratories/bessel_functions/hankel_toy.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy laboratories/bessel_functions/hankel_toy.cpp -lquadmath

$(TEST_BIN_DIR)/hankel_toy128: laboratories/bessel_functions/hankel_toy128.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy128 laboratories/bessel_functions/hankel_toy128.cpp -lquadmath

$(TEST_BIN_DIR)/hankel_toy_new: laboratories/bessel_functions/hankel_toy_new.cpp
	$(CXX17) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy_new laboratories/bessel_functions/hankel_toy_new.cpp -lquadmath


${CHECK_DIR}/check_airy_ai: ${CHECK_DIR}/check_airy_ai.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_ai ${CHECK_DIR}/check_airy_ai.cc -lquadmath

${CHECK_DIR}/check_airy_bi: ${CHECK_DIR}/check_airy_bi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_airy_bi ${CHECK_DIR}/check_airy_bi.cc -lquadmath

${CHECK_DIR}/check_assoc_laguerre: ${CHECK_DIR}/check_assoc_laguerre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_laguerre ${CHECK_DIR}/check_assoc_laguerre.cc -lquadmath

${CHECK_DIR}/check_assoc_legendre: ${CHECK_DIR}/check_assoc_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_assoc_legendre ${CHECK_DIR}/check_assoc_legendre.cc -lquadmath

${CHECK_DIR}/check_bernoulli: ${CHECK_DIR}/check_bernoulli.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_bernoulli ${CHECK_DIR}/check_bernoulli.cc -lquadmath

${CHECK_DIR}/check_beta: ${CHECK_DIR}/check_beta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_beta ${CHECK_DIR}/check_beta.cc -lquadmath

${CHECK_DIR}/check_binomial: ${CHECK_DIR}/check_binomial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_binomial ${CHECK_DIR}/check_binomial.cc -lquadmath

${CHECK_DIR}/check_chebyshev_t: ${CHECK_DIR}/check_chebyshev_t.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_t ${CHECK_DIR}/check_chebyshev_t.cc -lquadmath

${CHECK_DIR}/check_chebyshev_u: ${CHECK_DIR}/check_chebyshev_u.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_u ${CHECK_DIR}/check_chebyshev_u.cc -lquadmath

${CHECK_DIR}/check_chebyshev_v: ${CHECK_DIR}/check_chebyshev_v.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_v ${CHECK_DIR}/check_chebyshev_v.cc -lquadmath

${CHECK_DIR}/check_chebyshev_w: ${CHECK_DIR}/check_chebyshev_w.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chebyshev_w ${CHECK_DIR}/check_chebyshev_w.cc -lquadmath

${CHECK_DIR}/check_chi: ${CHECK_DIR}/check_chi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_chi ${CHECK_DIR}/check_chi.cc -lquadmath

${CHECK_DIR}/check_clausen_cl: ${CHECK_DIR}/check_clausen_cl.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_clausen_cl ${CHECK_DIR}/check_clausen_cl.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_1: ${CHECK_DIR}/check_comp_ellint_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_1 ${CHECK_DIR}/check_comp_ellint_1.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_2: ${CHECK_DIR}/check_comp_ellint_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_2 ${CHECK_DIR}/check_comp_ellint_2.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_3: ${CHECK_DIR}/check_comp_ellint_3.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_3 ${CHECK_DIR}/check_comp_ellint_3.cc -lquadmath

${CHECK_DIR}/check_comp_ellint_d: ${CHECK_DIR}/check_comp_ellint_d.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_comp_ellint_d ${CHECK_DIR}/check_comp_ellint_d.cc -lquadmath

${CHECK_DIR}/check_conf_hyperg: ${CHECK_DIR}/check_conf_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg ${CHECK_DIR}/check_conf_hyperg.cc -lquadmath

${CHECK_DIR}/check_conf_hyperg_lim: ${CHECK_DIR}/check_conf_hyperg_lim.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_conf_hyperg_lim ${CHECK_DIR}/check_conf_hyperg_lim.cc -lquadmath

${CHECK_DIR}/check_coshint: ${CHECK_DIR}/check_coshint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_coshint ${CHECK_DIR}/check_coshint.cc -lquadmath

${CHECK_DIR}/check_cosint: ${CHECK_DIR}/check_cosint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cosint ${CHECK_DIR}/check_cosint.cc -lquadmath

${CHECK_DIR}/check_cos_pi: ${CHECK_DIR}/check_cos_pi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cos_pi ${CHECK_DIR}/check_cos_pi.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_i: ${CHECK_DIR}/check_cyl_bessel_i.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_i ${CHECK_DIR}/check_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_j: ${CHECK_DIR}/check_cyl_bessel_j.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_j ${CHECK_DIR}/check_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/check_cyl_bessel_k: ${CHECK_DIR}/check_cyl_bessel_k.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_bessel_k ${CHECK_DIR}/check_cyl_bessel_k.cc -lquadmath

${CHECK_DIR}/check_cyl_hankel_1: ${CHECK_DIR}/check_cyl_hankel_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_1 ${CHECK_DIR}/check_cyl_hankel_1.cc -lquadmath

${CHECK_DIR}/check_cyl_hankel_2: ${CHECK_DIR}/check_cyl_hankel_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_hankel_2 ${CHECK_DIR}/check_cyl_hankel_2.cc -lquadmath

${CHECK_DIR}/check_cyl_neumann: ${CHECK_DIR}/check_cyl_neumann.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_cyl_neumann ${CHECK_DIR}/check_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_dawson: ${CHECK_DIR}/check_dawson.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dawson ${CHECK_DIR}/check_dawson.cc -lquadmath -lquadmath

${CHECK_DIR}/check_debye: ${CHECK_DIR}/check_debye.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_debye ${CHECK_DIR}/check_debye.cc -lquadmath -lquadmath

${CHECK_DIR}/check_digamma: ${CHECK_DIR}/check_digamma.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_digamma ${CHECK_DIR}/check_digamma.cc -lquadmath

${CHECK_DIR}/check_dilog: ${CHECK_DIR}/check_dilog.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dilog ${CHECK_DIR}/check_dilog.cc -lquadmath -lquadmath

${CHECK_DIR}/check_dirichlet_beta: ${CHECK_DIR}/check_dirichlet_beta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_beta ${CHECK_DIR}/check_dirichlet_beta.cc -lquadmath

${CHECK_DIR}/check_dirichlet_eta: ${CHECK_DIR}/check_dirichlet_eta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_eta ${CHECK_DIR}/check_dirichlet_eta.cc -lquadmath

${CHECK_DIR}/check_dirichlet_lambda: ${CHECK_DIR}/check_dirichlet_lambda.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_dirichlet_lambda ${CHECK_DIR}/check_dirichlet_lambda.cc -lquadmath

${CHECK_DIR}/check_double_factorial: ${CHECK_DIR}/check_double_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_double_factorial ${CHECK_DIR}/check_double_factorial.cc -lquadmath

${CHECK_DIR}/check_ellint_1: ${CHECK_DIR}/check_ellint_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_1 ${CHECK_DIR}/check_ellint_1.cc -lquadmath

${CHECK_DIR}/check_ellint_2: ${CHECK_DIR}/check_ellint_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_2 ${CHECK_DIR}/check_ellint_2.cc -lquadmath

${CHECK_DIR}/check_ellint_3: ${CHECK_DIR}/check_ellint_3.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_3 ${CHECK_DIR}/check_ellint_3.cc -lquadmath

${CHECK_DIR}/check_ellint_d: ${CHECK_DIR}/check_ellint_d.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_d ${CHECK_DIR}/check_ellint_d.cc -lquadmath

${CHECK_DIR}/check_ellint_rc: ${CHECK_DIR}/check_ellint_rc.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rc ${CHECK_DIR}/check_ellint_rc.cc -lquadmath

${CHECK_DIR}/check_ellint_rd: ${CHECK_DIR}/check_ellint_rd.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rd ${CHECK_DIR}/check_ellint_rd.cc -lquadmath

${CHECK_DIR}/check_ellint_rf: ${CHECK_DIR}/check_ellint_rf.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rf ${CHECK_DIR}/check_ellint_rf.cc -lquadmath

${CHECK_DIR}/check_ellint_rg: ${CHECK_DIR}/check_ellint_rg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rg ${CHECK_DIR}/check_ellint_rg.cc -lquadmath

${CHECK_DIR}/check_ellint_rj: ${CHECK_DIR}/check_ellint_rj.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellint_rj ${CHECK_DIR}/check_ellint_rj.cc -lquadmath

${CHECK_DIR}/check_ellnome: ${CHECK_DIR}/check_ellnome.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ellnome ${CHECK_DIR}/check_ellnome.cc -lquadmath

${CHECK_DIR}/check_euler: ${CHECK_DIR}/check_euler.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_euler ${CHECK_DIR}/check_euler.cc -lquadmath

${CHECK_DIR}/check_eulerian_1: ${CHECK_DIR}/check_eulerian_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_eulerian_1 ${CHECK_DIR}/check_eulerian_1.cc -lquadmath

${CHECK_DIR}/check_eulerian_2: ${CHECK_DIR}/check_eulerian_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_eulerian_2 ${CHECK_DIR}/check_eulerian_2.cc -lquadmath

${CHECK_DIR}/check_expint: ${CHECK_DIR}/check_expint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint ${CHECK_DIR}/check_expint.cc -lquadmath

${CHECK_DIR}/check_expint_en: ${CHECK_DIR}/check_expint_en.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_expint_en ${CHECK_DIR}/check_expint_en.cc -lquadmath

${CHECK_DIR}/check_factorial: ${CHECK_DIR}/check_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_factorial ${CHECK_DIR}/check_factorial.cc -lquadmath

${CHECK_DIR}/check_falling_factorial: ${CHECK_DIR}/check_falling_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_falling_factorial ${CHECK_DIR}/check_falling_factorial.cc -lquadmath

${CHECK_DIR}/check_fresnel_c: ${CHECK_DIR}/check_fresnel_c.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_c ${CHECK_DIR}/check_fresnel_c.cc -lquadmath

${CHECK_DIR}/check_fresnel_s: ${CHECK_DIR}/check_fresnel_s.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_fresnel_s ${CHECK_DIR}/check_fresnel_s.cc -lquadmath

${CHECK_DIR}/check_gamma_p: ${CHECK_DIR}/check_gamma_p.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_p ${CHECK_DIR}/check_gamma_p.cc -lquadmath

${CHECK_DIR}/check_gamma_q: ${CHECK_DIR}/check_gamma_q.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_q ${CHECK_DIR}/check_gamma_q.cc -lquadmath

${CHECK_DIR}/check_gamma_reciprocal: ${CHECK_DIR}/check_gamma_reciprocal.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gamma_reciprocal ${CHECK_DIR}/check_gamma_reciprocal.cc -lquadmath

${CHECK_DIR}/check_gegenbauer: ${CHECK_DIR}/check_gegenbauer.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_gegenbauer ${CHECK_DIR}/check_gegenbauer.cc -lquadmath

${CHECK_DIR}/check_hermite: ${CHECK_DIR}/check_hermite.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hermite ${CHECK_DIR}/check_hermite.cc -lquadmath

${CHECK_DIR}/check_heuman_lambda: ${CHECK_DIR}/check_heuman_lambda.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_heuman_lambda ${CHECK_DIR}/check_heuman_lambda.cc -lquadmath

${CHECK_DIR}/check_hurwitz_zeta: ${CHECK_DIR}/check_hurwitz_zeta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hurwitz_zeta ${CHECK_DIR}/check_hurwitz_zeta.cc -lquadmath

${CHECK_DIR}/check_hyperg: ${CHECK_DIR}/check_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_hyperg ${CHECK_DIR}/check_hyperg.cc -lquadmath

${CHECK_DIR}/check_ibeta: ${CHECK_DIR}/check_ibeta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ibeta ${CHECK_DIR}/check_ibeta.cc -lquadmath

${CHECK_DIR}/check_ibetac: ${CHECK_DIR}/check_ibetac.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ibetac ${CHECK_DIR}/check_ibetac.cc -lquadmath

${CHECK_DIR}/check_jacobi: ${CHECK_DIR}/check_jacobi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi ${CHECK_DIR}/check_jacobi.cc -lquadmath

${CHECK_DIR}/check_jacobi_cn: ${CHECK_DIR}/check_jacobi_cn.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_cn ${CHECK_DIR}/check_jacobi_cn.cc -lquadmath

${CHECK_DIR}/check_jacobi_dn: ${CHECK_DIR}/check_jacobi_dn.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_dn ${CHECK_DIR}/check_jacobi_dn.cc -lquadmath

${CHECK_DIR}/check_jacobi_sn: ${CHECK_DIR}/check_jacobi_sn.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_sn ${CHECK_DIR}/check_jacobi_sn.cc -lquadmath

${CHECK_DIR}/check_jacobi_zeta: ${CHECK_DIR}/check_jacobi_zeta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_jacobi_zeta ${CHECK_DIR}/check_jacobi_zeta.cc -lquadmath

${CHECK_DIR}/check_laguerre: ${CHECK_DIR}/check_laguerre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_laguerre ${CHECK_DIR}/check_laguerre.cc -lquadmath

${CHECK_DIR}/check_lbinomial: ${CHECK_DIR}/check_lbinomial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lbinomial ${CHECK_DIR}/check_lbinomial.cc -lquadmath

${CHECK_DIR}/check_ldouble_factorial: ${CHECK_DIR}/check_ldouble_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_ldouble_factorial ${CHECK_DIR}/check_ldouble_factorial.cc -lquadmath

${CHECK_DIR}/check_legendre: ${CHECK_DIR}/check_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre ${CHECK_DIR}/check_legendre.cc -lquadmath

${CHECK_DIR}/check_legendre_q: ${CHECK_DIR}/check_legendre_q.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_legendre_q ${CHECK_DIR}/check_legendre_q.cc -lquadmath

${CHECK_DIR}/check_lfactorial: ${CHECK_DIR}/check_lfactorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfactorial ${CHECK_DIR}/check_lfactorial.cc -lquadmath

${CHECK_DIR}/check_lfalling_factorial: ${CHECK_DIR}/check_lfalling_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lfalling_factorial ${CHECK_DIR}/check_lfalling_factorial.cc -lquadmath

${CHECK_DIR}/check_lgamma: ${CHECK_DIR}/check_lgamma.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lgamma ${CHECK_DIR}/check_lgamma.cc -lquadmath

${CHECK_DIR}/check_logistic_p: ${CHECK_DIR}/check_logistic_p.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_p ${CHECK_DIR}/check_logistic_p.cc -lquadmath

${CHECK_DIR}/check_logistic_pdf: ${CHECK_DIR}/check_logistic_pdf.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_logistic_pdf ${CHECK_DIR}/check_logistic_pdf.cc -lquadmath

${CHECK_DIR}/check_lognormal_p: ${CHECK_DIR}/check_lognormal_p.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_p ${CHECK_DIR}/check_lognormal_p.cc -lquadmath

${CHECK_DIR}/check_lognormal_pdf: ${CHECK_DIR}/check_lognormal_pdf.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lognormal_pdf ${CHECK_DIR}/check_lognormal_pdf.cc -lquadmath

${CHECK_DIR}/check_lrising_factorial: ${CHECK_DIR}/check_lrising_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_lrising_factorial ${CHECK_DIR}/check_lrising_factorial.cc -lquadmath

${CHECK_DIR}/check_normal_p: ${CHECK_DIR}/check_normal_p.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_p ${CHECK_DIR}/check_normal_p.cc -lquadmath

${CHECK_DIR}/check_normal_pdf: ${CHECK_DIR}/check_normal_pdf.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_normal_pdf ${CHECK_DIR}/check_normal_pdf.cc -lquadmath

${CHECK_DIR}/check_owens_t: ${CHECK_DIR}/check_owens_t.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_owens_t ${CHECK_DIR}/check_owens_t.cc -lquadmath

${CHECK_DIR}/check_polygamma: ${CHECK_DIR}/check_polygamma.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_polygamma ${CHECK_DIR}/check_polygamma.cc -lquadmath

${CHECK_DIR}/check_radpoly: ${CHECK_DIR}/check_radpoly.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_radpoly ${CHECK_DIR}/check_radpoly.cc -lquadmath

${CHECK_DIR}/check_riemann_zeta: ${CHECK_DIR}/check_riemann_zeta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_riemann_zeta ${CHECK_DIR}/check_riemann_zeta.cc -lquadmath

${CHECK_DIR}/check_rising_factorial: ${CHECK_DIR}/check_rising_factorial.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_rising_factorial ${CHECK_DIR}/check_rising_factorial.cc -lquadmath

${CHECK_DIR}/check_shi: ${CHECK_DIR}/check_shi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_shi ${CHECK_DIR}/check_shi.cc -lquadmath

${CHECK_DIR}/check_sinc: ${CHECK_DIR}/check_sinc.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc ${CHECK_DIR}/check_sinc.cc -lquadmath

${CHECK_DIR}/check_sinc_pi: ${CHECK_DIR}/check_sinc_pi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinc_pi ${CHECK_DIR}/check_sinc_pi.cc -lquadmath

${CHECK_DIR}/check_sinhint: ${CHECK_DIR}/check_sinhint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinhint ${CHECK_DIR}/check_sinhint.cc -lquadmath

${CHECK_DIR}/check_sinint: ${CHECK_DIR}/check_sinint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sinint ${CHECK_DIR}/check_sinint.cc -lquadmath

${CHECK_DIR}/check_sin_pi: ${CHECK_DIR}/check_sin_pi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sin_pi ${CHECK_DIR}/check_sin_pi.cc -lquadmath

${CHECK_DIR}/check_sph_bessel: ${CHECK_DIR}/check_sph_bessel.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel ${CHECK_DIR}/check_sph_bessel.cc -lquadmath

${CHECK_DIR}/check_sph_bessel_i: ${CHECK_DIR}/check_sph_bessel_i.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_i ${CHECK_DIR}/check_sph_bessel_i.cc -lquadmath

${CHECK_DIR}/check_sph_bessel_k: ${CHECK_DIR}/check_sph_bessel_k.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_bessel_k ${CHECK_DIR}/check_sph_bessel_k.cc -lquadmath

${CHECK_DIR}/check_sph_hankel_1: ${CHECK_DIR}/check_sph_hankel_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_1 ${CHECK_DIR}/check_sph_hankel_1.cc -lquadmath

${CHECK_DIR}/check_sph_hankel_2: ${CHECK_DIR}/check_sph_hankel_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_hankel_2 ${CHECK_DIR}/check_sph_hankel_2.cc -lquadmath

${CHECK_DIR}/check_sph_harmonic: ${CHECK_DIR}/check_sph_harmonic.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_harmonic ${CHECK_DIR}/check_sph_harmonic.cc -lquadmath

${CHECK_DIR}/check_sph_legendre: ${CHECK_DIR}/check_sph_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_legendre ${CHECK_DIR}/check_sph_legendre.cc -lquadmath

${CHECK_DIR}/check_sph_neumann: ${CHECK_DIR}/check_sph_neumann.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_sph_neumann ${CHECK_DIR}/check_sph_neumann.cc -lquadmath

${CHECK_DIR}/check_stirling_1: ${CHECK_DIR}/check_stirling_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_stirling_1 ${CHECK_DIR}/check_stirling_1.cc -lquadmath

${CHECK_DIR}/check_stirling_2: ${CHECK_DIR}/check_stirling_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_stirling_2 ${CHECK_DIR}/check_stirling_2.cc -lquadmath

${CHECK_DIR}/check_tgamma_lower: ${CHECK_DIR}/check_tgamma_lower.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma_lower ${CHECK_DIR}/check_tgamma_lower.cc -lquadmath

${CHECK_DIR}/check_tgamma: ${CHECK_DIR}/check_tgamma.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tgamma ${CHECK_DIR}/check_tgamma.cc -lquadmath

${CHECK_DIR}/check_theta_1: ${CHECK_DIR}/check_theta_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_1 ${CHECK_DIR}/check_theta_1.cc -lquadmath

${CHECK_DIR}/check_theta_2: ${CHECK_DIR}/check_theta_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_2 ${CHECK_DIR}/check_theta_2.cc -lquadmath

${CHECK_DIR}/check_theta_3: ${CHECK_DIR}/check_theta_3.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_3 ${CHECK_DIR}/check_theta_3.cc -lquadmath

${CHECK_DIR}/check_theta_4: ${CHECK_DIR}/check_theta_4.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_theta_4 ${CHECK_DIR}/check_theta_4.cc -lquadmath

${CHECK_DIR}/check_zernike: ${CHECK_DIR}/check_zernike.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_zernike ${CHECK_DIR}/check_zernike.cc -lquadmath

${CHECK_DIR}/complex_ellint_rc: ${CHECK_DIR}/complex_ellint_rc.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rc ${CHECK_DIR}/complex_ellint_rc.cc -lquadmath

${CHECK_DIR}/complex_ellint_rd: ${CHECK_DIR}/complex_ellint_rd.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rd ${CHECK_DIR}/complex_ellint_rd.cc -lquadmath

${CHECK_DIR}/complex_ellint_rf: ${CHECK_DIR}/complex_ellint_rf.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rf ${CHECK_DIR}/complex_ellint_rf.cc -lquadmath

${CHECK_DIR}/complex_ellint_rg: ${CHECK_DIR}/complex_ellint_rg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rg ${CHECK_DIR}/complex_ellint_rg.cc -lquadmath

${CHECK_DIR}/complex_ellint_rj: ${CHECK_DIR}/complex_ellint_rj.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_ellint_rj ${CHECK_DIR}/complex_ellint_rj.cc -lquadmath

${CHECK_DIR}/complex_airy_ai: ${CHECK_DIR}/complex_airy_ai.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_airy_ai ${CHECK_DIR}/complex_airy_ai.cc -lquadmath

${CHECK_DIR}/complex_airy_bi: ${CHECK_DIR}/complex_airy_bi.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/complex_airy_bi ${CHECK_DIR}/complex_airy_bi.cc -lquadmath

${CHECK_DIR}/deathmatch_comp_ellint: ${CHECK_DIR}/deathmatch_comp_ellint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_comp_ellint ${CHECK_DIR}/deathmatch_comp_ellint.cc -lquadmath

${CHECK_DIR}/deathmatch_conf_hyperg: ${CHECK_DIR}/deathmatch_conf_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_conf_hyperg ${CHECK_DIR}/deathmatch_conf_hyperg.cc -lquadmath

${CHECK_DIR}/deathmatch_conf_hyperg_lim: ${CHECK_DIR}/deathmatch_conf_hyperg_lim.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_conf_hyperg_lim ${CHECK_DIR}/deathmatch_conf_hyperg_lim.cc -lquadmath

${CHECK_DIR}/deathmatch_hyperg: ${CHECK_DIR}/deathmatch_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/deathmatch_hyperg ${CHECK_DIR}/deathmatch_hyperg.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_hankel_1: ${CHECK_DIR}/pr56216_cyl_hankel_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_1 ${CHECK_DIR}/pr56216_cyl_hankel_1.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_hankel_2: ${CHECK_DIR}/pr56216_cyl_hankel_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_hankel_2 ${CHECK_DIR}/pr56216_cyl_hankel_2.cc -lquadmath

${CHECK_DIR}/pr56216_cyl_bessel_i: ${CHECK_DIR}/pr56216_cyl_bessel_i.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr56216_cyl_bessel_i ${CHECK_DIR}/pr56216_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/pr68397: ${CHECK_DIR}/pr68397.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/pr68397 ${CHECK_DIR}/pr68397.cc -lquadmath

${CHECK_DIR}/origin_cyl_bessel_j: ${CHECK_DIR}/origin_cyl_bessel_j.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_bessel_j ${CHECK_DIR}/origin_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/origin_cyl_neumann: ${CHECK_DIR}/origin_cyl_neumann.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/origin_cyl_neumann ${CHECK_DIR}/origin_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_tr1_assoc_laguerre: ${CHECK_DIR}/check_tr1_assoc_laguerre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_assoc_laguerre ${CHECK_DIR}/check_tr1_assoc_laguerre.cc -lquadmath

${CHECK_DIR}/check_tr1_assoc_legendre: ${CHECK_DIR}/check_tr1_assoc_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_assoc_legendre ${CHECK_DIR}/check_tr1_assoc_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_beta: ${CHECK_DIR}/check_tr1_beta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_beta ${CHECK_DIR}/check_tr1_beta.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_1: ${CHECK_DIR}/check_tr1_comp_ellint_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_1 ${CHECK_DIR}/check_tr1_comp_ellint_1.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_2: ${CHECK_DIR}/check_tr1_comp_ellint_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_2 ${CHECK_DIR}/check_tr1_comp_ellint_2.cc -lquadmath

${CHECK_DIR}/check_tr1_comp_ellint_3: ${CHECK_DIR}/check_tr1_comp_ellint_3.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_comp_ellint_3 ${CHECK_DIR}/check_tr1_comp_ellint_3.cc -lquadmath

${CHECK_DIR}/check_tr1_conf_hyperg: ${CHECK_DIR}/check_tr1_conf_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_conf_hyperg ${CHECK_DIR}/check_tr1_conf_hyperg.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_i: ${CHECK_DIR}/check_tr1_cyl_bessel_i.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_i ${CHECK_DIR}/check_tr1_cyl_bessel_i.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_j: ${CHECK_DIR}/check_tr1_cyl_bessel_j.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_j ${CHECK_DIR}/check_tr1_cyl_bessel_j.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_bessel_k: ${CHECK_DIR}/check_tr1_cyl_bessel_k.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_bessel_k ${CHECK_DIR}/check_tr1_cyl_bessel_k.cc -lquadmath

${CHECK_DIR}/check_tr1_cyl_neumann: ${CHECK_DIR}/check_tr1_cyl_neumann.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_cyl_neumann ${CHECK_DIR}/check_tr1_cyl_neumann.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_1: ${CHECK_DIR}/check_tr1_ellint_1.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_1 ${CHECK_DIR}/check_tr1_ellint_1.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_2: ${CHECK_DIR}/check_tr1_ellint_2.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_2 ${CHECK_DIR}/check_tr1_ellint_2.cc -lquadmath

${CHECK_DIR}/check_tr1_ellint_3: ${CHECK_DIR}/check_tr1_ellint_3.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_ellint_3 ${CHECK_DIR}/check_tr1_ellint_3.cc -lquadmath

${CHECK_DIR}/check_tr1_expint: ${CHECK_DIR}/check_tr1_expint.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_expint ${CHECK_DIR}/check_tr1_expint.cc -lquadmath

${CHECK_DIR}/check_tr1_hermite: ${CHECK_DIR}/check_tr1_hermite.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_hermite ${CHECK_DIR}/check_tr1_hermite.cc -lquadmath

${CHECK_DIR}/check_tr1_hyperg: ${CHECK_DIR}/check_tr1_hyperg.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_hyperg ${CHECK_DIR}/check_tr1_hyperg.cc -lquadmath

${CHECK_DIR}/check_tr1_laguerre: ${CHECK_DIR}/check_tr1_laguerre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_laguerre ${CHECK_DIR}/check_tr1_laguerre.cc -lquadmath

${CHECK_DIR}/check_tr1_legendre: ${CHECK_DIR}/check_tr1_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_legendre ${CHECK_DIR}/check_tr1_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_riemann_zeta: ${CHECK_DIR}/check_tr1_riemann_zeta.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_riemann_zeta ${CHECK_DIR}/check_tr1_riemann_zeta.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_bessel: ${CHECK_DIR}/check_tr1_sph_bessel.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_bessel ${CHECK_DIR}/check_tr1_sph_bessel.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_legendre: ${CHECK_DIR}/check_tr1_sph_legendre.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_legendre ${CHECK_DIR}/check_tr1_sph_legendre.cc -lquadmath

${CHECK_DIR}/check_tr1_sph_neumann: ${CHECK_DIR}/check_tr1_sph_neumann.cc
	$(CXX14) -Iinclude -Ipolynomial/include -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o ${CHECK_DIR}/check_tr1_sph_neumann ${CHECK_DIR}/check_tr1_sph_neumann.cc -lquadmath


$(CHECK_DIR):
	if test ! -d $(CHECK_DIR); then \
	  mkdir $(CHECK_DIR); \
	fi


$(TEST_BIN_DIR):
	if test ! -d $(TEST_BIN_DIR); then \
	  mkdir $(TEST_BIN_DIR); \
	fi

$(TEST_OUT_DIR):
	if test ! -d $(TEST_OUT_DIR); then \
	  mkdir $(TEST_OUT_DIR); \
	fi

$(TEST_SF_OUT_DIR):
	if test ! -d $(TEST_SF_OUT_DIR); then \
	  mkdir $(TEST_SF_OUT_DIR); \
	fi

$(DIFF_SF_OUT_DIR):
	if test ! -d $(DIFF_SF_OUT_DIR); then \
	  mkdir $(DIFF_SF_OUT_DIR); \
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
	rm -f test/gsl_*_[fdlq].txt
	rm -f test/std_*_[fdlq].txt
	rm -f test/tr1_*_[fdlq].txt

diffclean:
	rm -f diff/diff_*_[fdlq].txt
