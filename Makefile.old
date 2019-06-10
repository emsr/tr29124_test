
# -Wconversion

CXX_VER = -std=gnu++17
CXX_INST_DIR = $(HOME)/bin$(SUFFIX)
ifeq ("$(wildcard $(CXX_INST_DIR)/bin/g++)","")
  CXX_VER = -std=gnu++2a
  CXX_INST_DIR = $(HOME)/bin
  ifeq ("$(wildcard $(CXX_INST_DIR)/bin/g++)","")
    CXX_VER = -std=gnu++17
    ifeq ($(wildcard "/mingw64"),"")
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
#CXXMAX = $(CXX17)
CXXMAX = $(CXX_INST_DIR)/bin/g++ $(CXX_VER) -g -Wall -Wextra -Wno-psabi
CXX_LIB_DIR = $(CXX_INST_DIR)/lib64
CXX_TEST_INC_DIR = libstdc++_support

INC_DIR = include/bits
INCLUDES = -Icxx_math_const/include -Iinclude -Icxx_fp_utils/include -Icxx_summation/include -Ipolynomial/include -Iquadrature/include
MPREAL_INCLUDES = -I$(HOME)/mpreal

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
       $(TEST_BIN_DIR)/plot_airy \
       $(TEST_BIN_DIR)/plot_bessel \
       $(TEST_BIN_DIR)/plot_gamma \
       $(TEST_BIN_DIR)/test_airy_roots \
       $(TEST_BIN_DIR)/test_anger_weber \
       $(TEST_BIN_DIR)/test_appell_f1 \
       $(TEST_BIN_DIR)/test_arith_geom_mean \
       $(TEST_BIN_DIR)/test_assoc_laguerre \
       $(TEST_BIN_DIR)/test_assoc_legendre \
       $(TEST_BIN_DIR)/test_assoc_legendre_q \
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
       $(TEST_BIN_DIR)/test_differentiation \
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
       $(TEST_BIN_DIR)/test_gegenbauer_neg_parm \
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
       $(TEST_BIN_DIR)/test_jacobi_neg_parm \
       $(TEST_BIN_DIR)/test_jacobi_neg_roots \
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
       $(TEST_BIN_DIR)/test_legendre_q \
       $(TEST_BIN_DIR)/test_legendre_ellint \
       $(TEST_BIN_DIR)/test_lentz_continued_fraction \
       $(TEST_BIN_DIR)/test_lerch \
       $(TEST_BIN_DIR)/test_limits \
       $(TEST_BIN_DIR)/test_little_airy \
       $(TEST_BIN_DIR)/test_lobatto \
       $(TEST_BIN_DIR)/test_log \
       $(TEST_BIN_DIR)/test_logsumexp \
       $(TEST_BIN_DIR)/test_lommel \
       $(TEST_BIN_DIR)/test_marcum_q \
       $(TEST_BIN_DIR)/test_math_h \
       $(TEST_BIN_DIR)/test_maxint \
       $(TEST_BIN_DIR)/test_meixner \
       $(TEST_BIN_DIR)/test_meixner_pollaczek \
       $(TEST_BIN_DIR)/test_mittag_leffler \
       $(TEST_BIN_DIR)/test_mod2pi \
       $(TEST_BIN_DIR)/test_mod_bessel_asymp \
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
       $(TEST_BIN_DIR)/test_ulp \
       $(TEST_BIN_DIR)/test_weierstrass_ellint \
       $(TEST_BIN_DIR)/test_wilson \
       $(TEST_BIN_DIR)/test_wright_omega \
       $(TEST_BIN_DIR)/test_zeta_trig \
       $(TEST_BIN_DIR)/run_coulfg \
       $(TEST_BIN_DIR)/RUN_COULFG \
       $(TEST_BIN_DIR)/build_bernoulli_2n_table \
       $(TEST_BIN_DIR)/build_zeta_trig_tables \
       $(TEST_BIN_DIR)/build_zeta_deriv_table \
       $(TEST_BIN_DIR)/build_etam1_table \
       $(TEST_BIN_DIR)/build_zetam1_table \
       $(TEST_BIN_DIR)/build_nfact_zetanp1 \
       $(TEST_BIN_DIR)/build_zetahalfm1_table \
       $(TEST_BIN_DIR)/build_gamma_lanczos \
       $(TEST_BIN_DIR)/build_gamma_spouge \
       $(TEST_BIN_DIR)/build_gamma_recip \
       $(TEST_BIN_DIR)/build_inv_erf_coefs \
       $(TEST_BIN_DIR)/build_sqrt_table \
       $(TEST_BIN_DIR)/build_sincos_tables \
       $(TEST_BIN_DIR)/build_atan_table \
       $(TEST_BIN_DIR)/build_cordic \
       $(TEST_BIN_DIR)/build_log_table


CHECKS = $(CHECK_DIR)/check_airy_ai \
	 $(CHECK_DIR)/check_airy_bi \
	 $(CHECK_DIR)/check_assoc_laguerre \
	 $(CHECK_DIR)/check_assoc_legendre \
	 $(CHECK_DIR)/check_bell \
	 $(CHECK_DIR)/check_beta \
	 $(CHECK_DIR)/check_bernoulli \
	 $(CHECK_DIR)/check_binomial \
	 $(CHECK_DIR)/check_chebyshev_t \
	 $(CHECK_DIR)/check_chebyshev_u \
	 $(CHECK_DIR)/check_chebyshev_v \
	 $(CHECK_DIR)/check_chebyshev_w \
	 $(CHECK_DIR)/check_chi \
	 $(CHECK_DIR)/check_clausen_cl \
	 $(CHECK_DIR)/check_comp_ellint_1 \
	 $(CHECK_DIR)/check_comp_ellint_2 \
	 $(CHECK_DIR)/check_comp_ellint_3 \
	 $(CHECK_DIR)/check_comp_ellint_d \
	 $(CHECK_DIR)/check_conf_hyperg \
	 $(CHECK_DIR)/check_conf_hyperg_lim \
	 $(CHECK_DIR)/check_coshint \
	 $(CHECK_DIR)/check_cosint \
	 $(CHECK_DIR)/check_cos_pi \
	 $(CHECK_DIR)/check_cyl_bessel_i \
	 $(CHECK_DIR)/check_cyl_bessel_i_scaled \
	 $(CHECK_DIR)/check_cyl_bessel_j \
	 $(CHECK_DIR)/check_cyl_bessel_k \
	 $(CHECK_DIR)/check_cyl_bessel_k_scaled \
	 $(CHECK_DIR)/check_cyl_hankel_1 \
	 $(CHECK_DIR)/check_cyl_hankel_2 \
	 $(CHECK_DIR)/check_cyl_neumann \
	 $(CHECK_DIR)/check_dawson \
	 $(CHECK_DIR)/check_debye \
	 $(CHECK_DIR)/check_digamma \
	 $(CHECK_DIR)/check_dilog \
	 $(CHECK_DIR)/check_dirichlet_beta \
	 $(CHECK_DIR)/check_dirichlet_eta \
	 $(CHECK_DIR)/check_dirichlet_lambda \
	 $(CHECK_DIR)/check_double_factorial \
	 $(CHECK_DIR)/check_ellint_1 \
	 $(CHECK_DIR)/check_ellint_2 \
	 $(CHECK_DIR)/check_ellint_3 \
	 $(CHECK_DIR)/check_ellint_d \
	 $(CHECK_DIR)/check_ellint_rc \
	 $(CHECK_DIR)/check_ellint_rd \
	 $(CHECK_DIR)/check_ellint_rf \
	 $(CHECK_DIR)/check_ellint_rg \
	 $(CHECK_DIR)/check_ellint_rj \
	 $(CHECK_DIR)/check_ellnome \
	 $(CHECK_DIR)/check_euler \
	 $(CHECK_DIR)/check_eulerian_1 \
	 $(CHECK_DIR)/check_eulerian_2 \
	 $(CHECK_DIR)/check_expint \
	 $(CHECK_DIR)/check_expint_en \
	 $(CHECK_DIR)/check_factorial \
	 $(CHECK_DIR)/check_falling_factorial \
	 $(CHECK_DIR)/check_fresnel_c \
	 $(CHECK_DIR)/check_fresnel_s \
	 $(CHECK_DIR)/check_gamma_p \
	 $(CHECK_DIR)/check_gamma_q \
	 $(CHECK_DIR)/check_gamma_reciprocal \
	 $(CHECK_DIR)/check_gegenbauer \
	 $(CHECK_DIR)/check_hermite \
	 $(CHECK_DIR)/check_heuman_lambda \
	 $(CHECK_DIR)/check_hurwitz_zeta \
	 $(CHECK_DIR)/check_hyperg \
	 $(CHECK_DIR)/check_ibeta \
	 $(CHECK_DIR)/check_ibetac \
	 $(CHECK_DIR)/check_jacobi \
	 $(CHECK_DIR)/check_jacobi_cn \
	 $(CHECK_DIR)/check_jacobi_dn \
	 $(CHECK_DIR)/check_jacobi_sn \
	 $(CHECK_DIR)/check_jacobi_zeta \
	 $(CHECK_DIR)/check_laguerre \
	 $(CHECK_DIR)/check_lah \
	 $(CHECK_DIR)/check_lbinomial \
	 $(CHECK_DIR)/check_ldouble_factorial \
	 $(CHECK_DIR)/check_legendre \
	 $(CHECK_DIR)/check_legendre \
	 $(CHECK_DIR)/check_legendre_q \
	 $(CHECK_DIR)/check_lfactorial \
	 $(CHECK_DIR)/check_lfalling_factorial \
	 $(CHECK_DIR)/check_lgamma \
	 $(CHECK_DIR)/check_logistic_p \
	 $(CHECK_DIR)/check_logistic_pdf \
	 $(CHECK_DIR)/check_lrising_factorial \
	 $(CHECK_DIR)/check_lognormal_p \
	 $(CHECK_DIR)/check_lognormal_pdf \
	 $(CHECK_DIR)/check_normal_p \
	 $(CHECK_DIR)/check_normal_pdf \
	 $(CHECK_DIR)/check_owens_t \
	 $(CHECK_DIR)/check_polygamma \
	 $(CHECK_DIR)/check_radpoly \
	 $(CHECK_DIR)/check_riemann_zeta \
	 $(CHECK_DIR)/check_rising_factorial \
	 $(CHECK_DIR)/check_shi \
	 $(CHECK_DIR)/check_sinc \
	 $(CHECK_DIR)/check_sinc_pi \
	 $(CHECK_DIR)/check_sinhint \
	 $(CHECK_DIR)/check_sinint \
	 $(CHECK_DIR)/check_sin_pi \
	 $(CHECK_DIR)/check_sph_bessel \
	 $(CHECK_DIR)/check_sph_bessel_i \
	 $(CHECK_DIR)/check_sph_bessel_k \
	 $(CHECK_DIR)/check_sph_hankel_1 \
	 $(CHECK_DIR)/check_sph_hankel_2 \
	 $(CHECK_DIR)/check_sph_harmonic \
	 $(CHECK_DIR)/check_sph_legendre \
	 $(CHECK_DIR)/check_sph_neumann \
	 $(CHECK_DIR)/check_stirling_1 \
	 $(CHECK_DIR)/check_stirling_2 \
	 $(CHECK_DIR)/check_tgamma_lower \
	 $(CHECK_DIR)/check_tgamma \
	 $(CHECK_DIR)/check_theta_1 \
	 $(CHECK_DIR)/check_theta_2 \
	 $(CHECK_DIR)/check_theta_3 \
	 $(CHECK_DIR)/check_theta_4 \
	 $(CHECK_DIR)/check_zernike \
	 $(CHECK_DIR)/complex_ellint_rc \
	 $(CHECK_DIR)/complex_ellint_rd \
	 $(CHECK_DIR)/complex_ellint_rf \
	 $(CHECK_DIR)/complex_ellint_rg \
	 $(CHECK_DIR)/complex_ellint_rj \
	 $(CHECK_DIR)/complex_airy_ai \
	 $(CHECK_DIR)/complex_airy_bi \
	 $(CHECK_DIR)/deathmatch_comp_ellint \
	 $(CHECK_DIR)/deathmatch_conf_hyperg \
	 $(CHECK_DIR)/deathmatch_conf_hyperg_lim \
	 $(CHECK_DIR)/deathmatch_hyperg \
	 $(CHECK_DIR)/pr56216_cyl_hankel_1 \
	 $(CHECK_DIR)/pr56216_cyl_hankel_2 \
	 $(CHECK_DIR)/pr56216_cyl_bessel_i \
	 $(CHECK_DIR)/pr68397 \
	 $(CHECK_DIR)/origin_cyl_bessel_j \
	 $(CHECK_DIR)/origin_cyl_neumann \
	 $(CHECK_DIR)/pr86655_assoc_legendre \
	 $(CHECK_DIR)/pr86655_sph_legendre

TR1_CHECKS =  \
	$(CHECK_DIR)/check_tr1_assoc_laguerre \
	$(CHECK_DIR)/check_tr1_assoc_legendre \
	$(CHECK_DIR)/check_tr1_beta \
	$(CHECK_DIR)/check_tr1_comp_ellint_1 \
	$(CHECK_DIR)/check_tr1_comp_ellint_2 \
	$(CHECK_DIR)/check_tr1_comp_ellint_3 \
	$(CHECK_DIR)/check_tr1_conf_hyperg \
	$(CHECK_DIR)/check_tr1_cyl_bessel_i \
	$(CHECK_DIR)/check_tr1_cyl_bessel_j \
	$(CHECK_DIR)/check_tr1_cyl_bessel_k \
	$(CHECK_DIR)/check_tr1_cyl_neumann \
	$(CHECK_DIR)/check_tr1_ellint_1 \
	$(CHECK_DIR)/check_tr1_ellint_2 \
	$(CHECK_DIR)/check_tr1_ellint_3 \
	$(CHECK_DIR)/check_tr1_expint \
	$(CHECK_DIR)/check_tr1_hermite \
	$(CHECK_DIR)/check_tr1_hyperg \
	$(CHECK_DIR)/check_tr1_laguerre \
	$(CHECK_DIR)/check_tr1_legendre \
	$(CHECK_DIR)/check_tr1_riemann_zeta \
	$(CHECK_DIR)/check_tr1_sph_bessel \
	$(CHECK_DIR)/check_tr1_sph_legendre \
	$(CHECK_DIR)/check_tr1_sph_neumann \
	$(CHECK_DIR)/pr86655_tr1_assoc_legendre \
	$(CHECK_DIR)/pr86655_tr1_sph_legendre


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


docs: Doxyfile include/bits/*
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
test: $(TEST_OUT_DIR) \
  run_airy_toy \
  run_hankel_toy \
  run_hankel_toy128 \
  run_hankel_toy_new \
  run_test_airy_roots \
  run_test_anger_weber \
  run_test_appell_f1 \
  run_test_assoc_laguerre \
  run_test_assoc_legendre \
  run_test_assoc_legendre_q \
  run_test_bernoulli \
  run_test_bessel \
  run_test_bessel_asymp \
  run_test_bessel_iter \
  run_test_beta \
  run_test_beta_inc \
  run_test_binet \
  run_test_binet_float \
  run_test_bose_einstein \
  run_test_charlier \
  run_test_chebyshev \
  run_test_chebyshev_trig \
  run_test_chebyshev_trig_pi \
  run_test_clausen \
  run_test_cmath \
  run_test_comp_ellint_1 \
  run_test_complex128 \
  run_test_complex_gamma \
  run_test_conf_hyperg \
  run_test_conf_hyperg_limit \
  run_test_const \
  run_test_continued_fraction \
  run_test_continuous_dual_hahn \
  run_test_continuous_hahn \
  run_test_cordic \
  run_test_coulomb \
  run_test_csint \
  run_test_cyl_hankel \
  run_test_dawson \
  run_test_debye \
  run_test_differentiation \
  run_test_digamma \
  run_test_dilog \
  run_test_dirichlet_eta \
  run_test_dual_hahn \
  run_test_erfc \
  run_test_experfc \
  run_test_expint \
  run_test_factorial \
  run_test_faddeeva \
  run_test_falling_factorial \
  run_test_fermi_dirac \
  run_test_fibonacci \
  run_test_float128 \
  run_test_fresnel \
  run_test_gamma \
  run_test_gamma_ratio \
  run_test_gamma_reciprocal \
  run_test_gegenbauer \
  run_test_gegenbauer_neg_parm \
  run_test_gudermannian \
  run_test_hahn \
  run_test_hankel \
  run_test_hankel_real_arg \
  run_test_hermite \
  run_test_heuman_lambda \
  run_test_hurwitz_zeta \
  run_test_hurwitz_zeta_new \
  run_test_hydrogen \
  run_test_hyperg \
  run_test_hypot \
  run_test_inv_erf \
  run_test_inv_gamma \
  run_test_inv_ibeta \
  run_test_inv_lgamma \
  run_test_jacobi \
  run_test_jacobi_neg_parm \
  run_test_jacobi_neg_roots \
  run_test_jacobi_ellint \
  run_test_jacobi_inv \
  run_test_jacobi_theta \
  run_test_jacobi_zeta \
  run_test_kelvin \
  run_test_krawtchouk \
  run_test_laguerre \
  run_test_lambert_w \
  run_test_legendre \
  run_test_legendre_q \
  run_test_legendre_ellint \
  run_test_lentz_continued_fraction \
  run_test_lerch \
  run_test_limits \
  run_test_little_airy \
  run_test_lobatto \
  run_test_log \
  run_test_logsumexp \
  run_test_lommel \
  run_test_marcum_q \
  run_test_math_h \
  run_test_maxint \
  run_test_meixner \
  run_test_meixner_pollaczek \
  run_test_mittag_leffler \
  run_test_mod2pi \
  run_test_mod_bessel_asymp \
  run_test_mpreal \
  run_test_notsospecfun \
  run_test_nric_bessel \
  run_test_numeric_limits \
  run_test_owens_t \
  run_test_parab_cyl \
  run_test_polygamma \
  run_test_polylog \
  run_test_power_mean \
  run_test_power_norm \
  run_test_pow_limits \
  run_test_racah \
  run_test_rational \
  run_test_recursion \
  run_test_reperiodized_hyper \
  run_test_reperiodized_trig \
  run_test_riemann_zeta \
  run_test_rising_factorial \
  run_test_root_search \
  run_test_sincos \
  run_test_sinus_cardinal \
  run_test_sph_bessel \
  run_test_sph_hankel \
  run_test_steed_continued_fraction \
  run_test_struve \
  run_test_struve_old \
  run_test_summation \
  run_test_theta \
  run_test_tr1_cmath \
  run_test_tricomi_u \
  run_test_trig \
  run_test_ulp \
  run_test_weierstrass_ellint \
  run_test_wilson \
  run_test_wright_omega \
  run_test_zeta_trig \
  run_run_coulfg \
  run_RUN_COULFG

plot: \
  run_airy_toy \
  run_airy_toy_old \
  run_plot_airy \
  run_plot_bessel \
  run_test_kelvin \
  run_test_struve \
  run_plot_gamma \
  run_test_beta \
  run_test_riemann_zeta

check: $(CHECK_DIR) $(CHECKS)
	echo "Beginning executions of checks..." > $(CHECK_DIR)/check_out.txt 2> $(CHECK_DIR)/check_err.txt
	echo "check_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bell" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_bell >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_bernoulli" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_bernoulli >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_binomial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_binomial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_chebyshev_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_u" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_chebyshev_u >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_v" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_chebyshev_v >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chebyshev_w" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_chebyshev_w >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_chi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_chi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_clausen_cl" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_clausen_cl >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_comp_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_comp_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_comp_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_comp_ellint_d" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_comp_ellint_d >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_conf_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_conf_hyperg_lim" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_conf_hyperg_lim >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_coshint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_coshint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cosint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cosint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cos_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cos_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_i_scaled" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_bessel_i_scaled >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_k" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_bessel_k_scaled" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_bessel_k_scaled >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dawson" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_dawson >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_debye" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_debye >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_digamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_digamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dilog" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_dilog >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_beta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_dirichlet_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_eta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_dirichlet_eta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_dirichlet_lambda" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_dirichlet_lambda >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_double_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_double_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_d" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_d >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_rc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rd" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_rd >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_rf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_rg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellint_rj" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellint_rj >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ellnome" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ellnome >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_euler" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_euler >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_eulerian_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_eulerian_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_eulerian_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_eulerian_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_expint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_expint_en" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_expint_en >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_falling_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_falling_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_c" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_fresnel_c >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_fresnel_s" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_fresnel_s >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_gamma_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_gamma_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gamma_reciprocal" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_gamma_reciprocal >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_gegenbauer" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_gegenbauer >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hermite" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_hermite >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_heuman_lambda" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_heuman_lambda >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hurwitz_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_hurwitz_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ibeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ibeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ibetac" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ibetac >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_jacobi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_cn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_jacobi_cn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_dn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_jacobi_dn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_sn" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_jacobi_sn >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_jacobi_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_jacobi_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_laguerre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lbinomial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lbinomial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lah" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lah >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_ldouble_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_ldouble_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_legendre_q" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_legendre_q >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lfactorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lfactorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_logistic_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_logistic_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_logistic_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lognormal_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lognormal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lognormal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lfalling_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lfalling_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_lrising_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_lrising_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_p" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_normal_p >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_normal_pdf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_normal_pdf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_owens_t" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_owens_t >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_polygamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_polygamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_radpoly" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_radpoly >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_riemann_zeta" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_rising_factorial" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_rising_factorial >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_shi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_shi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sinc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinc_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sinc_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinhint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sinhint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sinint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sinint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sin_pi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sin_pi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_bessel >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_bessel_k" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_harmonic" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_harmonic >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_sph_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_sph_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_stirling_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_stirling_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_stirling_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_stirling_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tgamma_lower" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tgamma_lower >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tgamma" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tgamma >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_theta_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_theta_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_3" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_theta_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_theta_4" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_theta_4 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_zernike" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_zernike >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rc" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_ellint_rc >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rd" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_ellint_rd >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rf" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_ellint_rf >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_ellint_rg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_ellint_rj" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_ellint_rj >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_airy_ai" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_airy_ai >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "complex_airy_bi" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/complex_airy_bi >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_comp_ellint" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/deathmatch_comp_ellint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_conf_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/deathmatch_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_conf_hyperg_lim" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/deathmatch_conf_hyperg_lim >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "deathmatch_hyperg" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/deathmatch_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_1" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr56216_cyl_hankel_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_hankel_2" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr56216_cyl_hankel_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr56216_cyl_bessel_i" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr56216_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr68397" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr68397 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_bessel_j" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/origin_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "origin_cyl_neumann" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/origin_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr86655_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr86655_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr86655_sph_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr86655_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


check_tr1: $(CHECK_DIR) $(TR1_CHECKS)
	echo "Beginning executions of tr1 checks..." > $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt
	echo "check_tr1_assoc_laguerre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_assoc_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_assoc_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_beta" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_beta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_1" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_comp_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_2" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_comp_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_comp_ellint_3" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_comp_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_conf_hyperg" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_conf_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_i" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_cyl_bessel_i >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_j" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_cyl_bessel_j >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_bessel_k" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_cyl_bessel_k >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_cyl_neumann" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_cyl_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_1" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_ellint_1 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_2" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_ellint_2 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_ellint_3" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_ellint_3 >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_expint" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_expint >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_hermite" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_hermite >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_hyperg" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_hyperg >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_laguerre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_laguerre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_riemann_zeta" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_riemann_zeta >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_bessel" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_sph_bessel >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_legendre" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "check_tr1_sph_neumann" >> $(CHECK_DIR)/check_tr1_out.txt 2> $(CHECK_DIR)/check_tr1_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/check_tr1_sph_neumann >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr86655_tr1_assoc_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr86655_tr1_assoc_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt
	echo "pr86655_tr1_sph_legendre" >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt && LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(CHECK_DIR)/pr86655_tr1_sph_legendre >> $(CHECK_DIR)/check_out.txt 2>> $(CHECK_DIR)/check_err.txt


libstdc++_support/testcase2: wrappers_debug libstdc++_support/testcase2.cpp libstdc++_support/testcase2.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXXMAX) $(INCLUDES) -o libstdc++_support/testcase2 libstdc++_support/testcase2.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXXMAX) -UTR1 $(INCLUDES) -o libstdc++_support/testcase libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

libstdc++_support/testcase_tr1: wrappers_debug libstdc++_support/testcase.cpp libstdc++_support/testcase.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXXMAX) -DTR1 $(INCLUDES) -o libstdc++_support/testcase_tr1 libstdc++_support/testcase.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

$(TEST_BIN_DIR)/plot_airy: laboratories/airy_functions/plot_airy.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/plot_airy laboratories/airy_functions/plot_airy.cpp -lquadmath

run_plot_airy: $(TEST_BIN_DIR)/plot_airy
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/plot_airy laboratories/plot_data > $(TEST_OUT_DIR)/plot_airy.txt

$(TEST_BIN_DIR)/plot_bessel: laboratories/bessel_functions/plot_bessel.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/plot_bessel laboratories/bessel_functions/plot_bessel.cpp -lquadmath

run_plot_bessel: $(TEST_BIN_DIR)/plot_bessel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/plot_bessel laboratories/plot_data > $(TEST_OUT_DIR)/plot_data.txt

$(TEST_BIN_DIR)/plot_gamma: laboratories/gamma_functions/plot_gamma.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/plot_gamma laboratories/gamma_functions/plot_gamma.cpp -lquadmath

run_plot_gamma: $(TEST_BIN_DIR)/plot_gamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/plot_gamma laboratories/plot_data > $(TEST_OUT_DIR)/plot_data.txt

$(TEST_BIN_DIR)/mpfrcalc: multiprecision/mpfr_gexpr.c
	$(GCC) -Iinclude -o $(TEST_BIN_DIR)/mpfrcalc multiprecision/mpfr_gexpr.c -lmpfr -lgmp -lm

$(TEST_BIN_DIR)/test_special_function: wrappers_debug laboratories/test_special_function.cpp laboratories/test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXXMAX) $(INCLUDES) -Ilaboratories -Iwrappers -o $(TEST_BIN_DIR)/test_special_function laboratories/test_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

$(TEST_BIN_DIR)/diff_special_function: wrappers_debug laboratories/diff_special_function.cpp laboratories/test_func.tcc $(INC_DIR)/*.h $(INC_DIR)/sf_*.tcc
	$(CXXMAX) $(INCLUDES) -Ilaboratories -Iwrappers -o $(TEST_BIN_DIR)/diff_special_function laboratories/diff_special_function.cpp -Wl,-rpath,$(CXX_LIB_DIR) -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost -lwrap_burkhardt -lgfortran -lwrap_cephes -lwrap_lerchphi -lwrap_faddeeva

$(TEST_BIN_DIR)/test_Faddeeva: $(FAD_DIR)/Faddeeva.hh $(FAD_DIR)/Faddeeva.cc
	$(CXX14) -DTEST_FADDEEVA -o $(TEST_BIN_DIR)/$(FAD_DIR)/test_Faddeeva $(FAD_DIR)/Faddeeva.cc -lquadmath

$(TEST_BIN_DIR)/test_airy_roots: laboratories/airy_functions/test_airy_roots.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_airy_roots laboratories/airy_functions/test_airy_roots.cpp -lquadmath

run_test_airy_roots: $(TEST_BIN_DIR)/test_airy_roots
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_airy_roots > $(TEST_OUT_DIR)/test_airy_roots.txt

$(TEST_BIN_DIR)/test_anger_weber: laboratories/bessel_functions/test_anger_weber.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_anger_weber laboratories/bessel_functions/test_anger_weber.cpp -lquadmath

run_test_anger_weber: $(TEST_BIN_DIR)/test_anger_weber
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_anger_weber > $(TEST_OUT_DIR)/test_anger_weber.txt

$(TEST_BIN_DIR)/test_appell_f1: laboratories/appell_functions/test_appell_f1.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_appell_f1 laboratories/appell_functions/test_appell_f1.cpp -lquadmath

run_test_appell_f1: $(TEST_BIN_DIR)/test_appell_f1
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_appell_f1 > $(TEST_OUT_DIR)/test_appell_f1.txt

$(TEST_BIN_DIR)/test_arith_geom_mean: laboratories/norm_functions/test_arith_geom_mean.cpp
	$(CXXMAX) -Iinclude -Icxx_fp_utils/include -Iwrappers -o $(TEST_BIN_DIR)/test_arith_geom_mean laboratories/norm_functions/test_arith_geom_mean.cpp -lquadmath

$(TEST_BIN_DIR)/test_assoc_laguerre: laboratories/orthogonal_polynomials/test_assoc_laguerre.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_assoc_laguerre laboratories/orthogonal_polynomials/test_assoc_laguerre.cpp -lquadmath

run_test_assoc_laguerre: $(TEST_BIN_DIR)/test_assoc_laguerre
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_assoc_laguerre > $(TEST_OUT_DIR)/test_assoc_laguerre.txt

$(TEST_BIN_DIR)/test_assoc_legendre: laboratories/orthogonal_polynomials/test_assoc_legendre.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_assoc_legendre laboratories/orthogonal_polynomials/test_assoc_legendre.cpp -lquadmath

run_test_assoc_legendre: $(TEST_BIN_DIR)/test_assoc_legendre
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_assoc_legendre > $(TEST_OUT_DIR)/test_assoc_legendre.txt

$(TEST_BIN_DIR)/test_assoc_legendre_q: laboratories/orthogonal_polynomials/test_assoc_legendre_q.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_assoc_legendre_q laboratories/orthogonal_polynomials/test_assoc_legendre_q.cpp -lquadmath

run_test_assoc_legendre_q: $(TEST_BIN_DIR)/test_assoc_legendre_q
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_assoc_legendre_q > $(TEST_OUT_DIR)/test_assoc_legendre_q.txt

$(TEST_BIN_DIR)/test_bernoulli: wrappers_debug laboratories/bernoulli_functions/test_bernoulli.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_bernoulli laboratories/bernoulli_functions/test_bernoulli.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_bernoulli: $(TEST_BIN_DIR)/test_bernoulli
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bernoulli > $(TEST_OUT_DIR)/test_bernoulli.txt

$(TEST_BIN_DIR)/test_bessel: laboratories/bessel_functions/test_bessel.cpp laboratories/bessel_functions/new_bessel.tcc
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel laboratories/bessel_functions/test_bessel.cpp -lquadmath

run_test_bessel: $(TEST_BIN_DIR)/test_bessel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel > $(TEST_OUT_DIR)/test_bessel.txt

$(TEST_BIN_DIR)/test_bessel_asymp: laboratories/bessel_functions/test_bessel_asymp.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel_asymp laboratories/bessel_functions/test_bessel_asymp.cpp -lquadmath

run_test_bessel_asymp: $(TEST_BIN_DIR)/test_bessel_asymp
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel_asymp > $(TEST_OUT_DIR)/test_bessel_asymp.txt

$(TEST_BIN_DIR)/test_bessel_iter: laboratories/bessel_functions/test_bessel_iter.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bessel_iter laboratories/bessel_functions/test_bessel_iter.cpp -lquadmath

run_test_bessel_iter: $(TEST_BIN_DIR)/test_bessel_iter
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bessel_iter > $(TEST_OUT_DIR)/test_bessel_iter.txt

$(TEST_BIN_DIR)/test_beta: wrappers_debug laboratories/beta_functions/test_beta.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_beta laboratories/beta_functions/test_beta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_beta: $(TEST_BIN_DIR)/test_beta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_beta laboratories/plot_data > $(TEST_OUT_DIR)/test_beta.txt

$(TEST_BIN_DIR)/test_beta_inc: laboratories/beta_functions/test_beta_inc.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_beta_inc laboratories/beta_functions/test_beta_inc.cpp -lquadmath

run_test_beta_inc: $(TEST_BIN_DIR)/test_beta_inc
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_beta_inc > $(TEST_OUT_DIR)/test_beta_inc.txt

$(TEST_BIN_DIR)/test_binet: laboratories/gamma_functions/test_binet.cpp
	$(CXXMAX) $(INCLUDES) -Icxx_rational/include -o $(TEST_BIN_DIR)/test_binet laboratories/gamma_functions/test_binet.cpp -lquadmath

run_test_binet: $(TEST_BIN_DIR)/test_binet
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_binet > $(TEST_OUT_DIR)/test_binet.txt

$(TEST_BIN_DIR)/test_binet_float: laboratories/gamma_functions/test_binet_float.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_binet_float laboratories/gamma_functions/test_binet_float.cpp -lquadmath

run_test_binet_float: $(TEST_BIN_DIR)/test_binet_float
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_binet_float > $(TEST_OUT_DIR)/test_binet_float.txt

$(TEST_BIN_DIR)/test_bose_einstein: laboratories/zeta_functions/test_bose_einstein.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_bose_einstein laboratories/zeta_functions/test_bose_einstein.cpp -lquadmath

run_test_bose_einstein: $(TEST_BIN_DIR)/test_bose_einstein
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_bose_einstein > $(TEST_OUT_DIR)/test_bose_einstein.txt

$(TEST_BIN_DIR)/test_charlier: laboratories/orthogonal_polynomials/test_charlier.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_charlier laboratories/orthogonal_polynomials/test_charlier.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_charlier: $(TEST_BIN_DIR)/test_charlier
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_charlier > $(TEST_OUT_DIR)/test_charlier.txt

$(TEST_BIN_DIR)/test_chebyshev: laboratories/orthogonal_polynomials/test_chebyshev.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev laboratories/orthogonal_polynomials/test_chebyshev.cpp -lquadmath

run_test_chebyshev: $(TEST_BIN_DIR)/test_chebyshev
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev > $(TEST_OUT_DIR)/test_chebyshev.txt

$(TEST_BIN_DIR)/test_chebyshev_trig: laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev_trig laboratories/orthogonal_polynomials/test_chebyshev_trig.cpp -lquadmath

run_test_chebyshev_trig: $(TEST_BIN_DIR)/test_chebyshev_trig
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev_trig > $(TEST_OUT_DIR)/test_chebyshev_trig.txt

$(TEST_BIN_DIR)/test_chebyshev_trig_pi: laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_chebyshev_trig_pi laboratories/orthogonal_polynomials/test_chebyshev_trig_pi.cpp -lquadmath

run_test_chebyshev_trig_pi: $(TEST_BIN_DIR)/test_chebyshev_trig_pi
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_chebyshev_trig_pi > $(TEST_OUT_DIR)/test_chebyshev_trig_pi.txt

$(TEST_BIN_DIR)/test_clausen: wrappers_debug laboratories/zeta_functions/test_clausen.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_clausen laboratories/zeta_functions/test_clausen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_clausen: $(TEST_BIN_DIR)/test_clausen
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_clausen > $(TEST_OUT_DIR)/test_clausen.txt

$(TEST_BIN_DIR)/test_cmath: test_std_maths/test_cmath.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_cmath test_std_maths/test_cmath.cpp -lquadmath

run_test_cmath: $(TEST_BIN_DIR)/test_cmath
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cmath > $(TEST_OUT_DIR)/test_cmath.txt

$(TEST_BIN_DIR)/test_comp_ellint_1: laboratories/elliptic_integrals/test_comp_ellint_1.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_comp_ellint_1 laboratories/elliptic_integrals/test_comp_ellint_1.cpp -lquadmath

run_test_comp_ellint_1: $(TEST_BIN_DIR)/test_comp_ellint_1
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_comp_ellint_1 > $(TEST_OUT_DIR)/test_comp_ellint_1.txt

$(TEST_BIN_DIR)/test_complex128: laboratories/complex_tools/test_complex128.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_complex128 laboratories/complex_tools/test_complex128.cpp -lquadmath

run_test_complex128: $(TEST_BIN_DIR)/test_complex128
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_complex128 > $(TEST_OUT_DIR)/test_complex128.txt

$(TEST_BIN_DIR)/test_complex_gamma: laboratories/gamma_functions/test_complex_gamma.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_complex_gamma laboratories/gamma_functions/test_complex_gamma.cpp -lquadmath

run_test_complex_gamma: $(TEST_BIN_DIR)/test_complex_gamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_complex_gamma > $(TEST_OUT_DIR)/test_complex_gamma.txt

$(TEST_BIN_DIR)/test_conf_hyperg: laboratories/hypergeometric_functions/test_conf_hyperg.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_conf_hyperg laboratories/hypergeometric_functions/test_conf_hyperg.cpp -lquadmath

run_test_conf_hyperg: $(TEST_BIN_DIR)/test_conf_hyperg
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_conf_hyperg > $(TEST_OUT_DIR)/test_conf_hyperg.txt

$(TEST_BIN_DIR)/test_conf_hyperg_limit: laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_conf_hyperg_limit laboratories/hypergeometric_functions/test_conf_hyperg_limit.cpp -lquadmath

run_test_conf_hyperg_limit: $(TEST_BIN_DIR)/test_conf_hyperg_limit
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_conf_hyperg_limit > $(TEST_OUT_DIR)/test_conf_hyperg_limit.txt

$(TEST_BIN_DIR)/test_const: laboratories/constants/test_const.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/test_const laboratories/constants/test_const.cpp -lquadmath -lmpfr -lgmp

run_test_const: $(TEST_BIN_DIR)/test_const
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_const > $(TEST_OUT_DIR)/test_const.txt

$(TEST_BIN_DIR)/test_continued_fraction: cxx_continued_fractions/test_continued_fraction.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continued_fraction cxx_continued_fractions/test_continued_fraction.cpp -lquadmath

run_test_continued_fraction: $(TEST_BIN_DIR)/test_continued_fraction
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continued_fraction > $(TEST_OUT_DIR)/test_continued_fraction.txt

$(TEST_BIN_DIR)/test_continuous_dual_hahn: laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continuous_dual_hahn laboratories/orthogonal_polynomials/test_continuous_dual_hahn.cpp -lquadmath

run_test_continuous_dual_hahn: $(TEST_BIN_DIR)/test_continuous_dual_hahn
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continuous_dual_hahn > $(TEST_OUT_DIR)/test_continuous_dual_hahn.txt

$(TEST_BIN_DIR)/test_continuous_hahn: laboratories/orthogonal_polynomials/test_continuous_hahn.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_continuous_hahn laboratories/orthogonal_polynomials/test_continuous_hahn.cpp -lquadmath

run_test_continuous_hahn: $(TEST_BIN_DIR)/test_continuous_hahn
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_continuous_hahn > $(TEST_OUT_DIR)/test_continuous_hahn.txt

$(TEST_BIN_DIR)/test_cordic: laboratories/elementary_functions/test_cordic.cpp
	$(CXXMAX) -Ilaboratories/elementary_functions $(INCLUDES) -o $(TEST_BIN_DIR)/test_cordic laboratories/elementary_functions/test_cordic.cpp -lquadmath

run_test_cordic: $(TEST_BIN_DIR)/test_cordic
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cordic > $(TEST_OUT_DIR)/test_cordic.txt

$(TEST_BIN_DIR)/test_coulomb: laboratories/coulomb_functions/test_coulomb.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_coulomb laboratories/coulomb_functions/test_coulomb.cpp -lquadmath

run_test_coulomb: $(TEST_BIN_DIR)/test_coulomb
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_coulomb > $(TEST_OUT_DIR)/test_coulomb.txt

$(TEST_BIN_DIR)/test_csint: laboratories/exponential_integrals/test_csint.cpp laboratories/exponential_integrals/csint.tcc
	$(CXXMAX) $(INCLUDES) -Ilaboratories/exponential_integrals -o $(TEST_BIN_DIR)/test_csint laboratories/exponential_integrals/test_csint.cpp -lquadmath

run_test_csint: $(TEST_BIN_DIR)/test_csint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_csint > $(TEST_OUT_DIR)/test_csint.txt

$(TEST_BIN_DIR)/test_cyl_hankel: wrappers_debug laboratories/bessel_functions/test_cyl_hankel.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_cyl_hankel laboratories/bessel_functions/test_cyl_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_cyl_hankel: $(TEST_BIN_DIR)/test_cyl_hankel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_cyl_hankel > $(TEST_OUT_DIR)/test_cyl_hankel.txt

$(TEST_BIN_DIR)/test_dawson: laboratories/error_functions/test_dawson.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dawson laboratories/error_functions/test_dawson.cpp -lquadmath

run_test_dawson: $(TEST_BIN_DIR)/test_dawson
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dawson > $(TEST_OUT_DIR)/test_dawson.txt

$(TEST_BIN_DIR)/test_debye: wrappers_debug laboratories/zeta_functions/test_debye.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_debye laboratories/zeta_functions/test_debye.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_debye: $(TEST_BIN_DIR)/test_debye
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_debye > $(TEST_OUT_DIR)/test_debye.txt

$(TEST_BIN_DIR)/test_differentiation: wrappers_debug cxx_differentiation/test_differentiation.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -Icxx_differentiation/include -o $(TEST_BIN_DIR)/test_differentiation cxx_differentiation/test_differentiation.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_differentiation: $(TEST_BIN_DIR)/test_differentiation
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_differentiation > $(TEST_OUT_DIR)/test_differentiation.txt

$(TEST_BIN_DIR)/test_digamma: wrappers_debug laboratories/gamma_functions/test_digamma.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_digamma laboratories/gamma_functions/test_digamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_digamma: $(TEST_BIN_DIR)/test_digamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_digamma > $(TEST_OUT_DIR)/test_digamma.txt

$(TEST_BIN_DIR)/test_dilog: laboratories/zeta_functions/test_dilog.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dilog laboratories/zeta_functions/test_dilog.cpp -lquadmath

run_test_dilog: $(TEST_BIN_DIR)/test_dilog
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dilog > $(TEST_OUT_DIR)/test_dilog.txt

$(TEST_BIN_DIR)/test_dirichlet_eta: laboratories/zeta_functions/test_dirichlet_eta.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dirichlet_eta laboratories/zeta_functions/test_dirichlet_eta.cpp -lquadmath

run_test_dirichlet_eta: $(TEST_BIN_DIR)/test_dirichlet_eta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dirichlet_eta > $(TEST_OUT_DIR)/test_dirichlet_eta.txt

$(TEST_BIN_DIR)/test_dual_hahn: laboratories/orthogonal_polynomials/test_dual_hahn.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_dual_hahn laboratories/orthogonal_polynomials/test_dual_hahn.cpp -lquadmath

run_test_dual_hahn: $(TEST_BIN_DIR)/test_dual_hahn
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_dual_hahn > $(TEST_OUT_DIR)/test_dual_hahn.txt

$(TEST_BIN_DIR)/test_erfc: wrappers_debug laboratories/error_functions/test_erfc.cpp
	$(CXXMAX) $(INCLUDES) -Icxx_continued_fractions/include -Iwrappers -o $(TEST_BIN_DIR)/test_erfc laboratories/error_functions/test_erfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_erfc: $(TEST_BIN_DIR)/test_erfc
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_erfc > $(TEST_OUT_DIR)/test_erfc.txt

$(TEST_BIN_DIR)/test_experfc: wrappers_debug laboratories/error_functions/test_experfc.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_experfc laboratories/error_functions/test_experfc.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost -lmpfr -lgmp

run_test_experfc: $(TEST_BIN_DIR)/test_experfc
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_experfc > $(TEST_OUT_DIR)/test_experfc.txt

$(TEST_BIN_DIR)/test_expint: wrappers_debug laboratories/exponential_integrals/test_expint.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_expint laboratories/exponential_integrals/test_expint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_expint: $(TEST_BIN_DIR)/test_expint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_expint > $(TEST_OUT_DIR)/test_expint.txt

$(TEST_BIN_DIR)/test_factorial: laboratories/gamma_functions/test_factorial.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_factorial laboratories/gamma_functions/test_factorial.cpp -lquadmath

run_test_factorial: $(TEST_BIN_DIR)/test_factorial
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_factorial > $(TEST_OUT_DIR)/test_factorial.txt

$(TEST_BIN_DIR)/test_faddeeva: wrappers_debug laboratories/error_functions/test_faddeeva.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_faddeeva laboratories/error_functions/test_faddeeva.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_faddeeva

run_test_faddeeva: $(TEST_BIN_DIR)/test_faddeeva
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_faddeeva > $(TEST_OUT_DIR)/test_faddeeva.txt

$(TEST_BIN_DIR)/test_falling_factorial: wrappers_debug laboratories/gamma_functions/test_falling_factorial.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_falling_factorial laboratories/gamma_functions/test_falling_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_falling_factorial: $(TEST_BIN_DIR)/test_falling_factorial
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_falling_factorial > $(TEST_OUT_DIR)/test_falling_factorial.txt

$(TEST_BIN_DIR)/test_fermi_dirac: wrappers_debug laboratories/zeta_functions/test_fermi_dirac.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_fermi_dirac laboratories/zeta_functions/test_fermi_dirac.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_fermi_dirac: $(TEST_BIN_DIR)/test_fermi_dirac
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fermi_dirac > $(TEST_OUT_DIR)/test_fermi_dirac.txt

$(TEST_BIN_DIR)/test_fibonacci: wrappers_debug laboratories/recursive_functions/test_fibonacci.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_fibonacci laboratories/recursive_functions/test_fibonacci.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_fibonacci: $(TEST_BIN_DIR)/test_fibonacci
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fibonacci > $(TEST_OUT_DIR)/test_fibonacci.txt

$(TEST_BIN_DIR)/test_float128: laboratories/floating_point_tools/test_float128.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_float128 laboratories/floating_point_tools/test_float128.cpp -lquadmath

run_test_float128: $(TEST_BIN_DIR)/test_float128
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_float128 > $(TEST_OUT_DIR)/test_float128.txt

$(TEST_BIN_DIR)/test_fresnel: wrappers_debug laboratories/error_functions/test_fresnel.cpp laboratories/error_functions/fresnel.tcc
	$(CXXMAX) $(INCLUDES) -Iwrappers -Ilaboratories/error_functions -o $(TEST_BIN_DIR)/test_fresnel laboratories/error_functions/test_fresnel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_fresnel: $(TEST_BIN_DIR)/test_fresnel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_fresnel > $(TEST_OUT_DIR)/test_fresnel.txt

$(TEST_BIN_DIR)/test_gamma: wrappers_debug laboratories/gamma_functions/test_gamma.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_gamma laboratories/gamma_functions/test_gamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_gamma: $(TEST_BIN_DIR)/test_gamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma > $(TEST_OUT_DIR)/test_gamma.txt

$(TEST_BIN_DIR)/test_gamma_ratio: wrappers_debug laboratories/gamma_functions/test_gamma_ratio.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_gamma_ratio laboratories/gamma_functions/test_gamma_ratio.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_gamma_ratio: $(TEST_BIN_DIR)/test_gamma_ratio
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma_ratio > $(TEST_OUT_DIR)/test_gamma_ratio.txt

$(TEST_BIN_DIR)/test_gamma_reciprocal: laboratories/gamma_functions/test_gamma_reciprocal.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_gamma_reciprocal laboratories/gamma_functions/test_gamma_reciprocal.cpp -lquadmath

run_test_gamma_reciprocal: $(TEST_BIN_DIR)/test_gamma_reciprocal
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gamma_reciprocal > $(TEST_OUT_DIR)/test_gamma_reciprocal.txt

$(TEST_BIN_DIR)/test_gegenbauer: laboratories/orthogonal_polynomials/test_gegenbauer.cpp
	$(CXXMAX) -Ilaboratories/orthogonal_polynomials $(INCLUDES) -o $(TEST_BIN_DIR)/test_gegenbauer laboratories/orthogonal_polynomials/test_gegenbauer.cpp -lquadmath

run_test_gegenbauer: $(TEST_BIN_DIR)/test_gegenbauer
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gegenbauer > $(TEST_OUT_DIR)/test_gegenbauer.txt

$(TEST_BIN_DIR)/test_gegenbauer_neg_parm: laboratories/orthogonal_polynomials/test_gegenbauer_neg_parm.cpp
	$(CXXMAX) $(INCLUDES) -Ilaboratories/orthogonal_polynomials -Iwrappers -o $(TEST_BIN_DIR)/test_gegenbauer_neg_parm laboratories/orthogonal_polynomials/test_gegenbauer_neg_parm.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_gegenbauer_neg_parm: $(TEST_BIN_DIR)/test_gegenbauer_neg_parm
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gegenbauer_neg_parm > $(TEST_OUT_DIR)/test_gegenbauer_neg_parm.txt

$(TEST_BIN_DIR)/test_gegenbauer_neg_roots: laboratories/orthogonal_polynomials/test_gegenbauer_neg_roots.cpp
	$(CXXMAX) $(INCLUDES) -Ilaboratories/orthogonal_polynomials -Ipolynomial/include -Iwrappers -o $(TEST_BIN_DIR)/test_gegenbauer_neg_roots laboratories/orthogonal_polynomials/test_gegenbauer_neg_roots.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_gegenbauer_neg_roots: $(TEST_BIN_DIR)/test_gegenbauer_neg_roots
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gegenbauer_neg_roots > $(TEST_OUT_DIR)/test_gegenbauer_neg_roots.txt

$(TEST_BIN_DIR)/test_gudermannian: laboratories/elementary_functions/test_gudermannian.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_gudermannian laboratories/elementary_functions/test_gudermannian.cpp -lquadmath

run_test_gudermannian: $(TEST_BIN_DIR)/test_gudermannian
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_gudermannian > $(TEST_OUT_DIR)/test_gudermannian.txt

$(TEST_BIN_DIR)/test_hahn: laboratories/orthogonal_polynomials/test_hahn.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hahn laboratories/orthogonal_polynomials/test_hahn.cpp -lquadmath

run_test_hahn: $(TEST_BIN_DIR)/test_hahn
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hahn > $(TEST_OUT_DIR)/test_hahn.txt

$(TEST_BIN_DIR)/test_hankel: laboratories/bessel_functions/test_hankel.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hankel laboratories/bessel_functions/test_hankel.cpp -lquadmath

run_test_hankel: $(TEST_BIN_DIR)/test_hankel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hankel > $(TEST_OUT_DIR)/test_hankel.txt

$(TEST_BIN_DIR)/test_hankel_real_arg: wrappers_debug laboratories/bessel_functions/test_hankel_real_arg.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hankel_real_arg laboratories/bessel_functions/test_hankel_real_arg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_hankel_real_arg: $(TEST_BIN_DIR)/test_hankel_real_arg
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hankel_real_arg > $(TEST_OUT_DIR)/test_hankel_real_arg.txt

$(TEST_BIN_DIR)/test_hermite: laboratories/orthogonal_polynomials/test_hermite.cpp laboratories/orthogonal_polynomials/new_hermite.tcc
	$(CXXMAX) $(INCLUDES) -Iquadrature/include -Icxx_continued_fractions/include -Ilaboratories/orthogonal_polynomials -o $(TEST_BIN_DIR)/test_hermite laboratories/orthogonal_polynomials/test_hermite.cpp -lquadmath

run_test_hermite: $(TEST_BIN_DIR)/test_hermite
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hermite > $(TEST_OUT_DIR)/test_hermite.txt

$(TEST_BIN_DIR)/test_heuman_lambda: wrappers_debug laboratories/elliptic_integrals/test_heuman_lambda.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_heuman_lambda laboratories/elliptic_integrals/test_heuman_lambda.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_heuman_lambda: $(TEST_BIN_DIR)/test_heuman_lambda
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_heuman_lambda > $(TEST_OUT_DIR)/test_heuman_lambda.txt

$(TEST_BIN_DIR)/test_hurwitz_zeta: laboratories/zeta_functions/test_hurwitz_zeta.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hurwitz_zeta laboratories/zeta_functions/test_hurwitz_zeta.cpp -lquadmath

run_test_hurwitz_zeta: $(TEST_BIN_DIR)/test_hurwitz_zeta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hurwitz_zeta > $(TEST_OUT_DIR)/test_hurwitz_zeta.txt

$(TEST_BIN_DIR)/test_hurwitz_zeta_new: laboratories/zeta_functions/test_hurwitz_zeta_new.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hurwitz_zeta_new laboratories/zeta_functions/test_hurwitz_zeta_new.cpp -lquadmath

run_test_hurwitz_zeta_new: $(TEST_BIN_DIR)/test_hurwitz_zeta_new
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hurwitz_zeta_new > $(TEST_OUT_DIR)/test_hurwitz_zeta_new.txt

$(TEST_BIN_DIR)/test_hydrogen: wrappers_debug laboratories/coulomb_functions/test_hydrogen.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hydrogen laboratories/coulomb_functions/test_hydrogen.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_hydrogen: $(TEST_BIN_DIR)/test_hydrogen
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hydrogen > $(TEST_OUT_DIR)/test_hydrogen.txt

$(TEST_BIN_DIR)/test_hyperg: wrappers_debug laboratories/hypergeometric_functions/test_hyperg.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_hyperg laboratories/hypergeometric_functions/test_hyperg.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_hyperg: $(TEST_BIN_DIR)/test_hyperg
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hyperg > $(TEST_OUT_DIR)/test_hyperg.txt

$(TEST_BIN_DIR)/test_hypot: laboratories/norm_functions/test_hypot.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_hypot laboratories/norm_functions/test_hypot.cpp -lquadmath

run_test_hypot: $(TEST_BIN_DIR)/test_hypot
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_hypot > $(TEST_OUT_DIR)/test_hypot.txt

$(TEST_BIN_DIR)/test_inv_erf: laboratories/error_functions/test_inv_erf.cpp
	$(CXXMAX) $(INCLUDES) -DSTANDALONE -o $(TEST_BIN_DIR)/test_inv_erf laboratories/error_functions/test_inv_erf.cpp -lquadmath

run_test_inv_erf: $(TEST_BIN_DIR)/test_inv_erf
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_erf > $(TEST_OUT_DIR)/test_inv_erf.txt

$(TEST_BIN_DIR)/test_inv_gamma: laboratories/gamma_functions/test_inv_gamma.cpp laboratories/error_functions/test_inv_erf.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_inv_gamma laboratories/gamma_functions/test_inv_gamma.cpp -lquadmath

run_test_inv_gamma: $(TEST_BIN_DIR)/test_inv_gamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_gamma > $(TEST_OUT_DIR)/test_inv_gamma.txt

$(TEST_BIN_DIR)/test_inv_ibeta: laboratories/beta_functions/test_inv_ibeta.cpp
	$(CXXMAX) $(INCLUDES) -Icxx_root_search/include -o $(TEST_BIN_DIR)/test_inv_ibeta laboratories/beta_functions/test_inv_ibeta.cpp -lquadmath

run_test_inv_ibeta: $(TEST_BIN_DIR)/test_inv_ibeta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_ibeta > $(TEST_OUT_DIR)/test_inv_ibeta.txt

$(TEST_BIN_DIR)/test_inv_lgamma: laboratories/gamma_functions/test_inv_lgamma.cpp
	$(CXXMAX) -Ilaboratories/gamma_functions $(INCLUDES) -o $(TEST_BIN_DIR)/test_inv_lgamma laboratories/gamma_functions/test_inv_lgamma.cpp -lquadmath

run_test_inv_lgamma: $(TEST_BIN_DIR)/test_inv_lgamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_inv_lgamma > $(TEST_OUT_DIR)/test_inv_lgamma.txt

$(TEST_BIN_DIR)/test_jacobi: laboratories/orthogonal_polynomials/test_jacobi.cpp
	$(CXXMAX) -Ilaboratories/orthogonal_polynomials $(INCLUDES) -o $(TEST_BIN_DIR)/test_jacobi laboratories/orthogonal_polynomials/test_jacobi.cpp -lquadmath

run_test_jacobi: $(TEST_BIN_DIR)/test_jacobi
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi > $(TEST_OUT_DIR)/test_jacobi.txt

$(TEST_BIN_DIR)/test_jacobi_neg_parm: laboratories/orthogonal_polynomials/test_jacobi_neg_parm.cpp
	$(CXXMAX) $(INCLUDES) -Ilaboratories/orthogonal_polynomials -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_neg_parm laboratories/orthogonal_polynomials/test_jacobi_neg_parm.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_jacobi_neg_parm: $(TEST_BIN_DIR)/test_jacobi_neg_parm
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_neg_parm > $(TEST_OUT_DIR)/test_jacobi_neg_parm.txt

$(TEST_BIN_DIR)/test_jacobi_neg_roots: laboratories/orthogonal_polynomials/test_jacobi_neg_roots.cpp
	$(CXXMAX) $(INCLUDES) -Ilaboratories/orthogonal_polynomials -Ipolynomial/include -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_neg_roots laboratories/orthogonal_polynomials/test_jacobi_neg_roots.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_jacobi_neg_roots: $(TEST_BIN_DIR)/test_jacobi_neg_roots
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_neg_roots > $(TEST_OUT_DIR)/test_jacobi_neg_roots.txt

$(TEST_BIN_DIR)/test_jacobi_ellint: laboratories/theta_functions/test_jacobi_ellint.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_ellint laboratories/theta_functions/test_jacobi_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

run_test_jacobi_ellint: $(TEST_BIN_DIR)/test_jacobi_ellint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_ellint > $(TEST_OUT_DIR)/test_jacobi_ellint.txt

$(TEST_BIN_DIR)/test_jacobi_inv: laboratories/theta_functions/test_jacobi_inv.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_jacobi_inv laboratories/theta_functions/test_jacobi_inv.cpp -lquadmath

run_test_jacobi_inv: $(TEST_BIN_DIR)/test_jacobi_inv
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_inv > $(TEST_OUT_DIR)/test_jacobi_inv.txt

$(TEST_BIN_DIR)/test_jacobi_theta: wrappers_debug laboratories/theta_functions/test_jacobi_theta.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_theta laboratories/theta_functions/test_jacobi_theta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_jacobi_theta: $(TEST_BIN_DIR)/test_jacobi_theta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_theta > $(TEST_OUT_DIR)/test_jacobi_theta.txt

$(TEST_BIN_DIR)/test_jacobi_zeta: wrappers_debug laboratories/elliptic_integrals/test_jacobi_zeta.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_jacobi_zeta laboratories/elliptic_integrals/test_jacobi_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_jacobi_zeta: $(TEST_BIN_DIR)/test_jacobi_zeta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_jacobi_zeta > $(TEST_OUT_DIR)/test_jacobi_zeta.txt

$(TEST_BIN_DIR)/test_kelvin: laboratories/bessel_functions/test_kelvin.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_kelvin laboratories/bessel_functions/test_kelvin.cpp -lquadmath

run_test_kelvin: $(TEST_BIN_DIR)/test_kelvin
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_kelvin laboratories/plot_data > $(TEST_OUT_DIR)/test_kelvin.txt

$(TEST_BIN_DIR)/test_krawtchouk: laboratories/orthogonal_polynomials/test_krawtchouk.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_krawtchouk laboratories/orthogonal_polynomials/test_krawtchouk.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_krawtchouk: $(TEST_BIN_DIR)/test_krawtchouk
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_krawtchouk > $(TEST_OUT_DIR)/test_krawtchouk.txt

$(TEST_BIN_DIR)/test_laguerre: laboratories/orthogonal_polynomials/test_laguerre.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_laguerre laboratories/orthogonal_polynomials/test_laguerre.cpp -lquadmath

run_test_laguerre: $(TEST_BIN_DIR)/test_laguerre
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_laguerre > $(TEST_OUT_DIR)/test_laguerre.txt

$(TEST_BIN_DIR)/test_lambert_w: laboratories/elementary_functions/test_lambert_w.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_lambert_w laboratories/elementary_functions/test_lambert_w.cpp -lquadmath

run_test_lambert_w: $(TEST_BIN_DIR)/test_lambert_w
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lambert_w > $(TEST_OUT_DIR)/test_lambert_w.txt

$(TEST_BIN_DIR)/test_large_order_bessel: laboratories/bessel_functions/test_large_order_bessel.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_large_order_bessel laboratories/bessel_functions/test_large_order_bessel.cpp -lquadmath

$(TEST_BIN_DIR)/test_legendre: laboratories/orthogonal_polynomials/test_legendre.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_legendre laboratories/orthogonal_polynomials/test_legendre.cpp -lquadmath

run_test_legendre: $(TEST_BIN_DIR)/test_legendre
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_legendre > $(TEST_OUT_DIR)/test_legendre.txt

$(TEST_BIN_DIR)/test_legendre_q: laboratories/orthogonal_polynomials/test_legendre_q.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_legendre_q laboratories/orthogonal_polynomials/test_legendre_q.cpp -lquadmath

run_test_legendre_q: $(TEST_BIN_DIR)/test_legendre_q
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_legendre_q > $(TEST_OUT_DIR)/test_legendre_q.txt

$(TEST_BIN_DIR)/test_legendre_ellint: laboratories/elliptic_integrals/test_legendre_ellint.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_legendre_ellint laboratories/elliptic_integrals/test_legendre_ellint.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_legendre_ellint: $(TEST_BIN_DIR)/test_legendre_ellint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_legendre_ellint > $(TEST_OUT_DIR)/test_legendre_ellint.txt

$(TEST_BIN_DIR)/test_lentz_continued_fraction: cxx_continued_fractions/test_lentz_continued_fraction.cpp
	$(CXXMAX) $(INCLUDES) -Icxx_continued_fractions/include -o $(TEST_BIN_DIR)/test_lentz_continued_fraction cxx_continued_fractions/test_lentz_continued_fraction.cpp -lquadmath

run_test_lentz_continued_fraction: $(TEST_BIN_DIR)/test_lentz_continued_fraction
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lentz_continued_fraction > $(TEST_OUT_DIR)/test_lentz_continued_fraction.txt

$(TEST_BIN_DIR)/test_lerch: laboratories/zeta_functions/test_lerch.cpp
	$(CXXMAX) $(INCLUDES) -I. -o $(TEST_BIN_DIR)/test_lerch laboratories/zeta_functions/test_lerch.cpp 3rdparty/lerchphi/Source/lerchphi.cpp -lquadmath

run_test_lerch: $(TEST_BIN_DIR)/test_lerch
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lerch > $(TEST_OUT_DIR)/test_lerch.txt

$(TEST_BIN_DIR)/test_limits: laboratories/floating_point_tools/test_limits.cpp
	$(CXXMAX) -Iinclude -o $(TEST_BIN_DIR)/test_limits laboratories/floating_point_tools/test_limits.cpp -lquadmath

run_test_limits: $(TEST_BIN_DIR)/test_limits
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_limits > $(TEST_OUT_DIR)/test_limits.txt

$(TEST_BIN_DIR)/test_little_airy: wrappers_debug laboratories/airy_functions/test_little_airy.cpp
	$(CXXMAX) -Ilaboratories/airy_functions $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_little_airy laboratories/airy_functions/test_little_airy.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_little_airy: $(TEST_BIN_DIR)/test_little_airy
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_little_airy > $(TEST_OUT_DIR)/test_little_airy.txt

$(TEST_BIN_DIR)/test_lobatto: wrappers_debug laboratories/orthogonal_polynomials/test_lobatto.cpp
	$(CXXMAX) -Ilaboratories/orthogonal_polynomials $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_lobatto laboratories/orthogonal_polynomials/test_lobatto.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_lobatto: $(TEST_BIN_DIR)/test_lobatto
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lobatto > $(TEST_OUT_DIR)/test_lobatto.txt

$(TEST_BIN_DIR)/test_log: laboratories/elementary_functions/test_log.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_log laboratories/elementary_functions/test_log.cpp -lquadmath

run_test_log: $(TEST_BIN_DIR)/test_log
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_log > $(TEST_OUT_DIR)/test_log.txt

$(TEST_BIN_DIR)/test_lommel: wrappers_debug laboratories/bessel_functions/test_lommel.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_lommel laboratories/bessel_functions/test_lommel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_lommel: $(TEST_BIN_DIR)/test_lommel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_lommel > $(TEST_OUT_DIR)/test_lommel.txt

$(TEST_BIN_DIR)/test_logsumexp: laboratories/norm_functions/test_logsumexp.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_logsumexp laboratories/norm_functions/test_logsumexp.cpp -lquadmath

run_test_logsumexp: $(TEST_BIN_DIR)/test_logsumexp
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_logsumexp > $(TEST_OUT_DIR)/test_logsumexp.txt

$(TEST_BIN_DIR)/test_marcum_q: laboratories/distribution_functions/test_marcum_q.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_marcum_q laboratories/distribution_functions/test_marcum_q.cpp -lquadmath

run_test_marcum_q: $(TEST_BIN_DIR)/test_marcum_q
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_marcum_q > $(TEST_OUT_DIR)/test_marcum_q.txt

$(TEST_BIN_DIR)/test_math_h: test_std_maths/test_math_h.cpp
	$(CXXMAX) $(INCLUDES) -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o $(TEST_BIN_DIR)/test_math_h test_std_maths/test_math_h.cpp -lquadmath

run_test_math_h: $(TEST_BIN_DIR)/test_math_h
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_math_h > $(TEST_OUT_DIR)/test_math_h.txt

$(TEST_BIN_DIR)/test_maxint: laboratories/floating_point_tools/test_maxint.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/test_maxint laboratories/floating_point_tools/test_maxint.cpp -lquadmath -lmpfr

run_test_maxint: $(TEST_BIN_DIR)/test_maxint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_maxint > $(TEST_OUT_DIR)/test_maxint.txt

$(TEST_BIN_DIR)/test_meixner: laboratories/orthogonal_polynomials/test_meixner.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_meixner laboratories/orthogonal_polynomials/test_meixner.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_meixner: $(TEST_BIN_DIR)/test_meixner
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_meixner > $(TEST_OUT_DIR)/test_meixner.txt

$(TEST_BIN_DIR)/test_meixner_pollaczek: laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_meixner_pollaczek laboratories/orthogonal_polynomials/test_meixner_pollaczek.cpp -lquadmath

run_test_meixner_pollaczek: $(TEST_BIN_DIR)/test_meixner_pollaczek
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_meixner_pollaczek > $(TEST_OUT_DIR)/test_meixner_pollaczek.txt

$(TEST_BIN_DIR)/test_mittag_leffler: laboratories/mittag_leffler_functions/test_mittag_leffler.cpp
	$(CXXMAX) $(INCLUDES) -Iquadrature/include -o $(TEST_BIN_DIR)/test_mittag_leffler laboratories/mittag_leffler_functions/test_mittag_leffler.cpp -lquadmath

run_test_mittag_leffler: $(TEST_BIN_DIR)/test_mittag_leffler
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mittag_leffler > $(TEST_OUT_DIR)/test_mittag_leffler.txt

$(TEST_BIN_DIR)/test_mod2pi: laboratories/floating_point_tools/test_mod2pi.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/test_mod2pi laboratories/floating_point_tools/test_mod2pi.cpp -lquadmath -lmpfr -lgmp

run_test_mod2pi: $(TEST_BIN_DIR)/test_mod2pi
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mod2pi > $(TEST_OUT_DIR)/test_mod2pi.txt

$(TEST_BIN_DIR)/test_mod_bessel_asymp: laboratories/bessel_functions/test_mod_bessel_asymp.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_mod_bessel_asymp laboratories/bessel_functions/test_mod_bessel_asymp.cpp -lquadmath

run_test_mod_bessel_asymp: $(TEST_BIN_DIR)/test_mod_bessel_asymp
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mod_bessel_asymp > $(TEST_OUT_DIR)/test_mod_bessel_asymp.txt

$(TEST_BIN_DIR)/test_mpreal: multiprecision/test_mpreal.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/test_mpreal multiprecision/test_mpreal.cpp -lquadmath -lmpfr -lgmp

run_test_mpreal: $(TEST_BIN_DIR)/test_mpreal
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_mpreal > $(TEST_OUT_DIR)/test_mpreal.txt

$(TEST_BIN_DIR)/test_notsospecfun: laboratories/elementary_functions/test_notsospecfun.cpp
	$(CXXMAX) -Iinclude -Icxx_fp_utils/include -o $(TEST_BIN_DIR)/test_notsospecfun laboratories/elementary_functions/test_notsospecfun.cpp -lquadmath

run_test_notsospecfun: $(TEST_BIN_DIR)/test_notsospecfun
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_notsospecfun > $(TEST_OUT_DIR)/test_notsospecfun.txt

$(TEST_BIN_DIR)/test_nric_bessel: laboratories/bessel_functions/test_nric_bessel.cpp laboratories/bessel_functions/nric_bessel.tcc
	$(CXXMAX) $(INCLUDES) -Ilaboratories/bessel_functions -o $(TEST_BIN_DIR)/test_nric_bessel laboratories/bessel_functions/test_nric_bessel.cpp -lquadmath

run_test_nric_bessel: $(TEST_BIN_DIR)/test_nric_bessel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_nric_bessel > $(TEST_OUT_DIR)/test_nric_bessel.txt

$(TEST_BIN_DIR)/test_numeric_limits: laboratories/floating_point_tools/test_numeric_limits.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/test_numeric_limits laboratories/floating_point_tools/test_numeric_limits.cpp -lquadmath -lmpfr -lgmp

run_test_numeric_limits: $(TEST_BIN_DIR)/test_numeric_limits
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_numeric_limits > $(TEST_OUT_DIR)/test_numeric_limits.txt

$(TEST_BIN_DIR)/test_owens_t: wrappers_debug laboratories/error_functions/test_owens_t.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_owens_t laboratories/error_functions/test_owens_t.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_owens_t: $(TEST_BIN_DIR)/test_owens_t
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_owens_t > $(TEST_OUT_DIR)/test_owens_t.txt

$(TEST_BIN_DIR)/test_parab_cyl: laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp
	$(CXXMAX) -Ilaboratories/parabolic_cylinder_functions $(INCLUDES) -o $(TEST_BIN_DIR)/test_parab_cyl laboratories/parabolic_cylinder_functions/test_parab_cyl.cpp -lquadmath

run_test_parab_cyl: $(TEST_BIN_DIR)/test_parab_cyl
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_parab_cyl > $(TEST_OUT_DIR)/test_parab_cyl.txt

$(TEST_BIN_DIR)/test_polygamma: laboratories/gamma_functions/test_polygamma.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -Icxx_continued_fractions/include -o $(TEST_BIN_DIR)/test_polygamma laboratories/gamma_functions/test_polygamma.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_polygamma: $(TEST_BIN_DIR)/test_polygamma
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_polygamma > $(TEST_OUT_DIR)/test_polygamma.txt

$(TEST_BIN_DIR)/test_polylog: wrappers_debug laboratories/zeta_functions/test_polylog.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_polylog laboratories/zeta_functions/test_polylog.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_cephes

run_test_polylog: $(TEST_BIN_DIR)/test_polylog
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_polylog > $(TEST_OUT_DIR)/test_polylog.txt

$(TEST_BIN_DIR)/test_power_mean: laboratories/norm_functions/test_power_mean.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_power_mean laboratories/norm_functions/test_power_mean.cpp -lquadmath

run_test_power_mean: $(TEST_BIN_DIR)/test_power_mean
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_power_mean > $(TEST_OUT_DIR)/test_power_mean.txt

$(TEST_BIN_DIR)/test_power_norm: laboratories/norm_functions/test_power_norm.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_power_norm laboratories/norm_functions/test_power_norm.cpp -lquadmath

run_test_power_norm: $(TEST_BIN_DIR)/test_power_norm
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_power_norm > $(TEST_OUT_DIR)/test_power_norm.txt

$(TEST_BIN_DIR)/test_pow_limits: test_std_maths/test_pow_limits.cpp
	$(CXXMAX) -DBIT -o $(TEST_BIN_DIR)/test_pow_limits test_std_maths/test_pow_limits.cpp -lquadmath

run_test_pow_limits: $(TEST_BIN_DIR)/test_pow_limits
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_pow_limits > $(TEST_OUT_DIR)/test_pow_limits.txt

$(TEST_BIN_DIR)/test_racah: laboratories/orthogonal_polynomials/test_racah.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_racah laboratories/orthogonal_polynomials/test_racah.cpp -lquadmath

run_test_racah: $(TEST_BIN_DIR)/test_racah
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_racah > $(TEST_OUT_DIR)/test_racah.txt

$(TEST_BIN_DIR)/test_rational: cxx_rational/test_rational.cpp
	$(CXXMAX) -Icxx_rational/include $(INCLUDES) -o $(TEST_BIN_DIR)/test_rational cxx_rational/test_rational.cpp -lquadmath

run_test_rational: $(TEST_BIN_DIR)/test_rational
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_rational > $(TEST_OUT_DIR)/test_rational.txt

$(TEST_BIN_DIR)/test_recursion: recursion/test_recursion.cpp
	$(CXXMAX) -Iinclude -o $(TEST_BIN_DIR)/test_recursion recursion/test_recursion.cpp -lquadmath

run_test_recursion: $(TEST_BIN_DIR)/test_recursion
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_recursion > $(TEST_OUT_DIR)/test_recursion.txt

$(TEST_BIN_DIR)/test_reperiodized_hyper: laboratories/elementary_functions/test_reperiodized_hyper.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_reperiodized_hyper laboratories/elementary_functions/test_reperiodized_hyper.cpp -lquadmath

run_test_reperiodized_hyper: $(TEST_BIN_DIR)/test_reperiodized_hyper
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_reperiodized_hyper > $(TEST_OUT_DIR)/test_reperiodized_hyper.txt

$(TEST_BIN_DIR)/test_reperiodized_trig: wrappers_debug laboratories/elementary_functions/test_reperiodized_trig.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_reperiodized_trig laboratories/elementary_functions/test_reperiodized_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_reperiodized_trig: $(TEST_BIN_DIR)/test_reperiodized_trig
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_reperiodized_trig > $(TEST_OUT_DIR)/test_reperiodized_trig.txt

$(TEST_BIN_DIR)/test_riemann_zeta: laboratories/zeta_functions/test_riemann_zeta.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_riemann_zeta laboratories/zeta_functions/test_riemann_zeta.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_riemann_zeta: $(TEST_BIN_DIR)/test_riemann_zeta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_riemann_zeta laboratories/plot_data > $(TEST_OUT_DIR)/test_riemann_zeta.txt

$(TEST_BIN_DIR)/test_rising_factorial: wrappers_debug laboratories/gamma_functions/test_rising_factorial.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_rising_factorial laboratories/gamma_functions/test_rising_factorial.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_rising_factorial: $(TEST_BIN_DIR)/test_rising_factorial
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_rising_factorial > $(TEST_OUT_DIR)/test_rising_factorial.txt

$(TEST_BIN_DIR)/test_root_search: cxx_root_search/test_root_search.cpp cxx_root_search/include/ext/*
	$(CXXMAX) $(INCLUDES) -Icxx_root_search/include -o $(TEST_BIN_DIR)/test_root_search cxx_root_search/test_root_search.cpp -lquadmath

run_test_root_search: $(TEST_BIN_DIR)/test_root_search
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_root_search > $(TEST_OUT_DIR)/test_root_search.txt

$(TEST_BIN_DIR)/test_sincos: laboratories/elementary_functions/test_sincos.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_sincos laboratories/elementary_functions/test_sincos.cpp -lquadmath

run_test_sincos: $(TEST_BIN_DIR)/test_sincos
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sincos > $(TEST_OUT_DIR)/test_sincos.txt

$(TEST_BIN_DIR)/test_sinus_cardinal: wrappers_debug laboratories/elementary_functions/test_sinus_cardinal.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_sinus_cardinal laboratories/elementary_functions/test_sinus_cardinal.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl -lwrap_boost

run_test_sinus_cardinal: $(TEST_BIN_DIR)/test_sinus_cardinal
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sinus_cardinal > $(TEST_OUT_DIR)/test_sinus_cardinal.txt

$(TEST_BIN_DIR)/test_sph_bessel: laboratories/bessel_functions/test_sph_bessel.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_sph_bessel laboratories/bessel_functions/test_sph_bessel.cpp -lquadmath

run_test_sph_bessel: $(TEST_BIN_DIR)/test_sph_bessel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sph_bessel > $(TEST_OUT_DIR)/test_sph_bessel.txt

$(TEST_BIN_DIR)/test_sph_hankel: wrappers_debug laboratories/bessel_functions/test_sph_hankel.cpp
	$(CXXMAX) $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_sph_hankel laboratories/bessel_functions/test_sph_hankel.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_boost

run_test_sph_hankel: $(TEST_BIN_DIR)/test_sph_hankel
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_sph_hankel > $(TEST_OUT_DIR)/test_sph_hankel.txt

$(TEST_BIN_DIR)/test_steed_continued_fraction: cxx_continued_fractions/test_steed_continued_fraction.cpp
	$(CXXMAX) $(INCLUDES) -Icxx_continued_fractions/include -o $(TEST_BIN_DIR)/test_steed_continued_fraction cxx_continued_fractions/test_steed_continued_fraction.cpp -lquadmath

run_test_steed_continued_fraction: $(TEST_BIN_DIR)/test_steed_continued_fraction
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_steed_continued_fraction > $(TEST_OUT_DIR)/test_steed_continued_fraction.txt

$(TEST_BIN_DIR)/test_struve: laboratories/bessel_functions/test_struve.cpp
	$(CXXMAX) -Ilaboratories/bessel_functions $(INCLUDES) -Iwrappers -o $(TEST_BIN_DIR)/test_struve laboratories/bessel_functions/test_struve.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_burkhardt -lgfortran

run_test_struve: $(TEST_BIN_DIR)/test_struve
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_struve laboratories/plot_data > $(TEST_OUT_DIR)/test_struve.txt

$(TEST_BIN_DIR)/test_struve_old: laboratories/bessel_functions/test_struve_old.cpp
	$(CXXMAX) -Ilaboratories/bessel_functions $(INCLUDES) -Ilaboratories/hypergeometric_functions -o $(TEST_BIN_DIR)/test_struve_old laboratories/bessel_functions/test_struve_old.cpp -lquadmath

run_test_struve_old: $(TEST_BIN_DIR)/test_struve_old
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_struve_old > $(TEST_OUT_DIR)/test_struve_old.txt

$(TEST_BIN_DIR)/test_summation: cxx_summation/test_summation.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_summation cxx_summation/test_summation.cpp -lquadmath

run_test_summation: $(TEST_BIN_DIR)/test_summation
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_summation > $(TEST_OUT_DIR)/test_summation.txt

$(TEST_BIN_DIR)/test_theta: laboratories/theta_functions/test_theta.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_theta laboratories/theta_functions/test_theta.cpp -lquadmath

run_test_theta: $(TEST_BIN_DIR)/test_theta
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_theta > $(TEST_OUT_DIR)/test_theta.txt

$(TEST_BIN_DIR)/test_tr1_cmath: test_std_maths/test_tr1_cmath.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_tr1_cmath test_std_maths/test_tr1_cmath.cpp -lquadmath

run_test_tr1_cmath: $(TEST_BIN_DIR)/test_tr1_cmath
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_tr1_cmath > $(TEST_OUT_DIR)/test_tr1_cmath.txt

$(TEST_BIN_DIR)/test_tricomi_u: laboratories/hypergeometric_functions/test_tricomi_u.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_tricomi_u laboratories/hypergeometric_functions/test_tricomi_u.cpp -lquadmath

run_test_tricomi_u: $(TEST_BIN_DIR)/test_tricomi_u
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_tricomi_u > $(TEST_OUT_DIR)/test_tricomi_u.txt

$(TEST_BIN_DIR)/test_trig: laboratories/elementary_functions/test_trig.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_trig laboratories/elementary_functions/test_trig.cpp -lquadmath

run_test_trig: $(TEST_BIN_DIR)/test_trig
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_trig > $(TEST_OUT_DIR)/test_trig.txt

$(TEST_BIN_DIR)/test_ulp: cxx_fp_utils/test_ulp.cpp
	$(CXXMAX) -Icxx_fp_utils/include -o $(TEST_BIN_DIR)/test_ulp cxx_fp_utils/test_ulp.cpp

run_test_ulp: $(TEST_BIN_DIR)/test_ulp
	$(TEST_BIN_DIR)/test_ulp > $(TEST_OUT_DIR)/test_ulp.txt

$(TEST_BIN_DIR)/test_weierstrass_ellint: laboratories/theta_functions/test_weierstrass_ellint.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_weierstrass_ellint laboratories/theta_functions/test_weierstrass_ellint.cpp -lquadmath

run_test_weierstrass_ellint: $(TEST_BIN_DIR)/test_weierstrass_ellint
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_weierstrass_ellint > $(TEST_OUT_DIR)/test_weierstrass_ellint.txt

$(TEST_BIN_DIR)/test_wilson: laboratories/orthogonal_polynomials/test_wilson.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_wilson laboratories/orthogonal_polynomials/test_wilson.cpp -lquadmath

run_test_wilson: $(TEST_BIN_DIR)/test_wilson
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_wilson > $(TEST_OUT_DIR)/test_wilson.txt

$(TEST_BIN_DIR)/test_wright_omega: laboratories/elementary_functions/test_wright_omega.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_wright_omega laboratories/elementary_functions/test_wright_omega.cpp -lquadmath

run_test_wright_omega: $(TEST_BIN_DIR)/test_wright_omega
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_wright_omega > $(TEST_OUT_DIR)/test_wright_omega.txt

$(TEST_BIN_DIR)/test_zeta_trig: laboratories/elementary_functions/test_zeta_trig.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/test_zeta_trig laboratories/elementary_functions/test_zeta_trig.cpp -lquadmath -L$(WRAP_DEBUG_DIR) -lwrap_gsl

run_test_zeta_trig: $(TEST_BIN_DIR)/test_zeta_trig
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/test_zeta_trig > $(TEST_OUT_DIR)/test_zeta_trig.txt


$(TEST_BIN_DIR)/run_coulfg: laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/run_coulfg laboratories/coulomb_functions/coulfg.cpp laboratories/coulomb_functions/run_coulfg.cpp -lquadmath

run_run_coulfg: $(TEST_BIN_DIR)/run_coulfg
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$(WRAP_DEBUG_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/run_coulfg > $(TEST_OUT_DIR)/run_coulfg.txt

$(TEST_BIN_DIR)/RUN_COULFG: laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR
	gfortran -o $(TEST_BIN_DIR)/RUN_COULFG laboratories/coulomb_functions/COULFG.FOR laboratories/coulomb_functions/RUN_COULFG.FOR

run_RUN_COULFG: ./RUN_COULFG
	./RUN_COULFG > $(TEST_OUT_DIR)/RUN_COULFG.TXT


$(TEST_BIN_DIR)/airy_toy: laboratories/airy_functions/airy_toy.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/airy_toy laboratories/airy_functions/airy_toy.cpp -lquadmath

run_airy_toy: $(TEST_BIN_DIR)/airy_toy
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/airy_toy laboratories/plot_data > $(TEST_OUT_DIR)/airy_toy.txt

$(TEST_BIN_DIR)/airy_toy_old: laboratories/airy_functions/airy_toy.cpp
	$(CXXMAX) $(INCLUDES) -DOLD -o $(TEST_BIN_DIR)/airy_toy_old laboratories/airy_functions/airy_toy.cpp -lquadmath

run_airy_toy_old: $(TEST_BIN_DIR)/airy_toy_old
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/airy_toy_old laboratories/plot_data > $(TEST_OUT_DIR)/airy_toy_old.txt

$(TEST_BIN_DIR)/hankel_toy: laboratories/bessel_functions/hankel_toy.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy laboratories/bessel_functions/hankel_toy.cpp -lquadmath

run_hankel_toy: $(TEST_BIN_DIR)/hankel_toy
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy > $(TEST_OUT_DIR)/hankel_toy.txt

$(TEST_BIN_DIR)/hankel_toy128: laboratories/bessel_functions/hankel_toy128.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy128 laboratories/bessel_functions/hankel_toy128.cpp -lquadmath

run_hankel_toy128: $(TEST_BIN_DIR)/hankel_toy128
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy128 > $(TEST_OUT_DIR)/hankel_toy128.txt

$(TEST_BIN_DIR)/hankel_toy_new: laboratories/bessel_functions/hankel_toy_new.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/hankel_toy_new laboratories/bessel_functions/hankel_toy_new.cpp -lquadmath

run_hankel_toy_new: $(TEST_BIN_DIR)/hankel_toy_new
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(TEST_BIN_DIR)/hankel_toy_new > $(TEST_OUT_DIR)/hankel_toy_new.txt


$(TEST_BIN_DIR)/build_bernoulli_2n_table: laboratories/bernoulli_functions/build_bernoulli_2n_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_bernoulli_2n_table laboratories/bernoulli_functions/build_bernoulli_2n_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_zeta_trig_tables: laboratories/zeta_functions/build_zeta_trig_tables.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_zeta_trig_tables laboratories/zeta_functions/build_zeta_trig_tables.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_zeta_deriv_table: laboratories/zeta_functions/build_zeta_deriv_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_zeta_deriv_table laboratories/zeta_functions/build_zeta_deriv_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_etam1_table: laboratories/zeta_functions/build_etam1_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_etam1_table laboratories/zeta_functions/build_etam1_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_zetam1_table: laboratories/zeta_functions/build_zetam1_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_zetam1_table laboratories/zeta_functions/build_zetam1_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_nfact_zetanp1: laboratories/zeta_functions/build_nfact_zetanp1.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_nfact_zetanp1 laboratories/zeta_functions/build_nfact_zetanp1.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_zetahalfm1_table: laboratories/zeta_functions/build_zetahalfm1_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_zetahalfm1_table laboratories/zeta_functions/build_zetahalfm1_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_gamma_lanczos: laboratories/gamma_functions/build_gamma_lanczos.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/build_gamma_lanczos laboratories/gamma_functions/build_gamma_lanczos.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_gamma_spouge: laboratories/gamma_functions/build_gamma_spouge.cpp
	$(CXXMAX) $(INCLUDES) -o $(TEST_BIN_DIR)/build_gamma_spouge laboratories/gamma_functions/build_gamma_spouge.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_gamma_recip: laboratories/gamma_functions/build_gamma_recip.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_gamma_recip laboratories/gamma_functions/build_gamma_recip.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_inv_erf_coefs: laboratories/error_functions/build_inv_erf_coefs.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_inv_erf_coefs laboratories/error_functions/build_inv_erf_coefs.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_sqrt_table: laboratories/elementary_functions/build_sqrt_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_sqrt_table laboratories/elementary_functions/build_sqrt_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_sincos_tables: laboratories/elementary_functions/build_sincos_tables.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_sincos_tables laboratories/elementary_functions/build_sincos_tables.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_atan_table: laboratories/elementary_functions/build_atan_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_atan_table laboratories/elementary_functions/build_atan_table.cpp -lquadmath -lmpfr

$(TEST_BIN_DIR)/build_cordic: laboratories/elementary_functions/build_cordic.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_cordic laboratories/elementary_functions/build_cordic.cpp -lquadmath -lmpfr -lgmp

$(TEST_BIN_DIR)/build_log_table: laboratories/elementary_functions/build_log_table.cpp
	$(CXXMAX) $(INCLUDES) $(MPREAL_INCLUDES) -o $(TEST_BIN_DIR)/build_log_table laboratories/elementary_functions/build_log_table.cpp -lquadmath -lmpfr


$(CHECK_DIR)/check_airy_ai: $(CHECK_DIR)/check_airy_ai.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_airy_ai $(CHECK_DIR)/check_airy_ai.cc -lquadmath

$(CHECK_DIR)/check_airy_bi: $(CHECK_DIR)/check_airy_bi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_airy_bi $(CHECK_DIR)/check_airy_bi.cc -lquadmath

$(CHECK_DIR)/check_assoc_laguerre: $(CHECK_DIR)/check_assoc_laguerre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_assoc_laguerre $(CHECK_DIR)/check_assoc_laguerre.cc -lquadmath

$(CHECK_DIR)/check_assoc_legendre: $(CHECK_DIR)/check_assoc_legendre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_assoc_legendre $(CHECK_DIR)/check_assoc_legendre.cc -lquadmath

$(CHECK_DIR)/check_bell: $(CHECK_DIR)/check_bell.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_bell $(CHECK_DIR)/check_bell.cc -lquadmath

$(CHECK_DIR)/check_bernoulli: $(CHECK_DIR)/check_bernoulli.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_bernoulli $(CHECK_DIR)/check_bernoulli.cc -lquadmath

$(CHECK_DIR)/check_beta: $(CHECK_DIR)/check_beta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_beta $(CHECK_DIR)/check_beta.cc -lquadmath

$(CHECK_DIR)/check_binomial: $(CHECK_DIR)/check_binomial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_binomial $(CHECK_DIR)/check_binomial.cc -lquadmath

$(CHECK_DIR)/check_chebyshev_t: $(CHECK_DIR)/check_chebyshev_t.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_chebyshev_t $(CHECK_DIR)/check_chebyshev_t.cc -lquadmath

$(CHECK_DIR)/check_chebyshev_u: $(CHECK_DIR)/check_chebyshev_u.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_chebyshev_u $(CHECK_DIR)/check_chebyshev_u.cc -lquadmath

$(CHECK_DIR)/check_chebyshev_v: $(CHECK_DIR)/check_chebyshev_v.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_chebyshev_v $(CHECK_DIR)/check_chebyshev_v.cc -lquadmath

$(CHECK_DIR)/check_chebyshev_w: $(CHECK_DIR)/check_chebyshev_w.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_chebyshev_w $(CHECK_DIR)/check_chebyshev_w.cc -lquadmath

$(CHECK_DIR)/check_chi: $(CHECK_DIR)/check_chi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_chi $(CHECK_DIR)/check_chi.cc -lquadmath

$(CHECK_DIR)/check_clausen_cl: $(CHECK_DIR)/check_clausen_cl.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_clausen_cl $(CHECK_DIR)/check_clausen_cl.cc -lquadmath

$(CHECK_DIR)/check_comp_ellint_1: $(CHECK_DIR)/check_comp_ellint_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_comp_ellint_1 $(CHECK_DIR)/check_comp_ellint_1.cc -lquadmath

$(CHECK_DIR)/check_comp_ellint_2: $(CHECK_DIR)/check_comp_ellint_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_comp_ellint_2 $(CHECK_DIR)/check_comp_ellint_2.cc -lquadmath

$(CHECK_DIR)/check_comp_ellint_3: $(CHECK_DIR)/check_comp_ellint_3.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_comp_ellint_3 $(CHECK_DIR)/check_comp_ellint_3.cc -lquadmath

$(CHECK_DIR)/check_comp_ellint_d: $(CHECK_DIR)/check_comp_ellint_d.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_comp_ellint_d $(CHECK_DIR)/check_comp_ellint_d.cc -lquadmath

$(CHECK_DIR)/check_conf_hyperg: $(CHECK_DIR)/check_conf_hyperg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_conf_hyperg $(CHECK_DIR)/check_conf_hyperg.cc -lquadmath

$(CHECK_DIR)/check_conf_hyperg_lim: $(CHECK_DIR)/check_conf_hyperg_lim.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_conf_hyperg_lim $(CHECK_DIR)/check_conf_hyperg_lim.cc -lquadmath

$(CHECK_DIR)/check_coshint: $(CHECK_DIR)/check_coshint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_coshint $(CHECK_DIR)/check_coshint.cc -lquadmath

$(CHECK_DIR)/check_cosint: $(CHECK_DIR)/check_cosint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cosint $(CHECK_DIR)/check_cosint.cc -lquadmath

$(CHECK_DIR)/check_cos_pi: $(CHECK_DIR)/check_cos_pi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cos_pi $(CHECK_DIR)/check_cos_pi.cc -lquadmath

$(CHECK_DIR)/check_cyl_bessel_i: $(CHECK_DIR)/check_cyl_bessel_i.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_bessel_i $(CHECK_DIR)/check_cyl_bessel_i.cc -lquadmath

$(CHECK_DIR)/check_cyl_bessel_i_scaled: $(CHECK_DIR)/check_cyl_bessel_i_scaled.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_bessel_i_scaled $(CHECK_DIR)/check_cyl_bessel_i_scaled.cc -lquadmath

$(CHECK_DIR)/check_cyl_bessel_j: $(CHECK_DIR)/check_cyl_bessel_j.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_bessel_j $(CHECK_DIR)/check_cyl_bessel_j.cc -lquadmath

$(CHECK_DIR)/check_cyl_bessel_k: $(CHECK_DIR)/check_cyl_bessel_k.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_bessel_k $(CHECK_DIR)/check_cyl_bessel_k.cc -lquadmath

$(CHECK_DIR)/check_cyl_bessel_k_scaled: $(CHECK_DIR)/check_cyl_bessel_k_scaled.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_bessel_k_scaled $(CHECK_DIR)/check_cyl_bessel_k_scaled.cc -lquadmath

$(CHECK_DIR)/check_cyl_hankel_1: $(CHECK_DIR)/check_cyl_hankel_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_hankel_1 $(CHECK_DIR)/check_cyl_hankel_1.cc -lquadmath

$(CHECK_DIR)/check_cyl_hankel_2: $(CHECK_DIR)/check_cyl_hankel_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_hankel_2 $(CHECK_DIR)/check_cyl_hankel_2.cc -lquadmath

$(CHECK_DIR)/check_cyl_neumann: $(CHECK_DIR)/check_cyl_neumann.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_cyl_neumann $(CHECK_DIR)/check_cyl_neumann.cc -lquadmath

$(CHECK_DIR)/check_dawson: $(CHECK_DIR)/check_dawson.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_dawson $(CHECK_DIR)/check_dawson.cc -lquadmath -lquadmath

$(CHECK_DIR)/check_debye: $(CHECK_DIR)/check_debye.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_debye $(CHECK_DIR)/check_debye.cc -lquadmath -lquadmath

$(CHECK_DIR)/check_digamma: $(CHECK_DIR)/check_digamma.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_digamma $(CHECK_DIR)/check_digamma.cc -lquadmath

$(CHECK_DIR)/check_dilog: $(CHECK_DIR)/check_dilog.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_dilog $(CHECK_DIR)/check_dilog.cc -lquadmath -lquadmath

$(CHECK_DIR)/check_dirichlet_beta: $(CHECK_DIR)/check_dirichlet_beta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_dirichlet_beta $(CHECK_DIR)/check_dirichlet_beta.cc -lquadmath

$(CHECK_DIR)/check_dirichlet_eta: $(CHECK_DIR)/check_dirichlet_eta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_dirichlet_eta $(CHECK_DIR)/check_dirichlet_eta.cc -lquadmath

$(CHECK_DIR)/check_dirichlet_lambda: $(CHECK_DIR)/check_dirichlet_lambda.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_dirichlet_lambda $(CHECK_DIR)/check_dirichlet_lambda.cc -lquadmath

$(CHECK_DIR)/check_double_factorial: $(CHECK_DIR)/check_double_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_double_factorial $(CHECK_DIR)/check_double_factorial.cc -lquadmath

$(CHECK_DIR)/check_ellint_1: $(CHECK_DIR)/check_ellint_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_1 $(CHECK_DIR)/check_ellint_1.cc -lquadmath

$(CHECK_DIR)/check_ellint_2: $(CHECK_DIR)/check_ellint_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_2 $(CHECK_DIR)/check_ellint_2.cc -lquadmath

$(CHECK_DIR)/check_ellint_3: $(CHECK_DIR)/check_ellint_3.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_3 $(CHECK_DIR)/check_ellint_3.cc -lquadmath

$(CHECK_DIR)/check_ellint_d: $(CHECK_DIR)/check_ellint_d.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_d $(CHECK_DIR)/check_ellint_d.cc -lquadmath

$(CHECK_DIR)/check_ellint_rc: $(CHECK_DIR)/check_ellint_rc.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_rc $(CHECK_DIR)/check_ellint_rc.cc -lquadmath

$(CHECK_DIR)/check_ellint_rd: $(CHECK_DIR)/check_ellint_rd.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_rd $(CHECK_DIR)/check_ellint_rd.cc -lquadmath

$(CHECK_DIR)/check_ellint_rf: $(CHECK_DIR)/check_ellint_rf.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_rf $(CHECK_DIR)/check_ellint_rf.cc -lquadmath

$(CHECK_DIR)/check_ellint_rg: $(CHECK_DIR)/check_ellint_rg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_rg $(CHECK_DIR)/check_ellint_rg.cc -lquadmath

$(CHECK_DIR)/check_ellint_rj: $(CHECK_DIR)/check_ellint_rj.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellint_rj $(CHECK_DIR)/check_ellint_rj.cc -lquadmath

$(CHECK_DIR)/check_ellnome: $(CHECK_DIR)/check_ellnome.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ellnome $(CHECK_DIR)/check_ellnome.cc -lquadmath

$(CHECK_DIR)/check_euler: $(CHECK_DIR)/check_euler.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_euler $(CHECK_DIR)/check_euler.cc -lquadmath

$(CHECK_DIR)/check_eulerian_1: $(CHECK_DIR)/check_eulerian_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_eulerian_1 $(CHECK_DIR)/check_eulerian_1.cc -lquadmath

$(CHECK_DIR)/check_eulerian_2: $(CHECK_DIR)/check_eulerian_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_eulerian_2 $(CHECK_DIR)/check_eulerian_2.cc -lquadmath

$(CHECK_DIR)/check_expint: $(CHECK_DIR)/check_expint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_expint $(CHECK_DIR)/check_expint.cc -lquadmath

$(CHECK_DIR)/check_expint_en: $(CHECK_DIR)/check_expint_en.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_expint_en $(CHECK_DIR)/check_expint_en.cc -lquadmath

$(CHECK_DIR)/check_factorial: $(CHECK_DIR)/check_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_factorial $(CHECK_DIR)/check_factorial.cc -lquadmath

$(CHECK_DIR)/check_falling_factorial: $(CHECK_DIR)/check_falling_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_falling_factorial $(CHECK_DIR)/check_falling_factorial.cc -lquadmath

$(CHECK_DIR)/check_fresnel_c: $(CHECK_DIR)/check_fresnel_c.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_fresnel_c $(CHECK_DIR)/check_fresnel_c.cc -lquadmath

$(CHECK_DIR)/check_fresnel_s: $(CHECK_DIR)/check_fresnel_s.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_fresnel_s $(CHECK_DIR)/check_fresnel_s.cc -lquadmath

$(CHECK_DIR)/check_gamma_p: $(CHECK_DIR)/check_gamma_p.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_gamma_p $(CHECK_DIR)/check_gamma_p.cc -lquadmath

$(CHECK_DIR)/check_gamma_q: $(CHECK_DIR)/check_gamma_q.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_gamma_q $(CHECK_DIR)/check_gamma_q.cc -lquadmath

$(CHECK_DIR)/check_gamma_reciprocal: $(CHECK_DIR)/check_gamma_reciprocal.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_gamma_reciprocal $(CHECK_DIR)/check_gamma_reciprocal.cc -lquadmath

$(CHECK_DIR)/check_gegenbauer: $(CHECK_DIR)/check_gegenbauer.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_gegenbauer $(CHECK_DIR)/check_gegenbauer.cc -lquadmath

$(CHECK_DIR)/check_hermite: $(CHECK_DIR)/check_hermite.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_hermite $(CHECK_DIR)/check_hermite.cc -lquadmath

$(CHECK_DIR)/check_heuman_lambda: $(CHECK_DIR)/check_heuman_lambda.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_heuman_lambda $(CHECK_DIR)/check_heuman_lambda.cc -lquadmath

$(CHECK_DIR)/check_hurwitz_zeta: $(CHECK_DIR)/check_hurwitz_zeta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_hurwitz_zeta $(CHECK_DIR)/check_hurwitz_zeta.cc -lquadmath

$(CHECK_DIR)/check_hyperg: $(CHECK_DIR)/check_hyperg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_hyperg $(CHECK_DIR)/check_hyperg.cc -lquadmath

$(CHECK_DIR)/check_ibeta: $(CHECK_DIR)/check_ibeta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ibeta $(CHECK_DIR)/check_ibeta.cc -lquadmath

$(CHECK_DIR)/check_ibetac: $(CHECK_DIR)/check_ibetac.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ibetac $(CHECK_DIR)/check_ibetac.cc -lquadmath

$(CHECK_DIR)/check_jacobi: $(CHECK_DIR)/check_jacobi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_jacobi $(CHECK_DIR)/check_jacobi.cc -lquadmath

$(CHECK_DIR)/check_jacobi_cn: $(CHECK_DIR)/check_jacobi_cn.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_jacobi_cn $(CHECK_DIR)/check_jacobi_cn.cc -lquadmath

$(CHECK_DIR)/check_jacobi_dn: $(CHECK_DIR)/check_jacobi_dn.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_jacobi_dn $(CHECK_DIR)/check_jacobi_dn.cc -lquadmath

$(CHECK_DIR)/check_jacobi_sn: $(CHECK_DIR)/check_jacobi_sn.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_jacobi_sn $(CHECK_DIR)/check_jacobi_sn.cc -lquadmath

$(CHECK_DIR)/check_jacobi_zeta: $(CHECK_DIR)/check_jacobi_zeta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_jacobi_zeta $(CHECK_DIR)/check_jacobi_zeta.cc -lquadmath

$(CHECK_DIR)/check_laguerre: $(CHECK_DIR)/check_laguerre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_laguerre $(CHECK_DIR)/check_laguerre.cc -lquadmath

$(CHECK_DIR)/check_lah: $(CHECK_DIR)/check_lah.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lah $(CHECK_DIR)/check_lah.cc -lquadmath

$(CHECK_DIR)/check_lbinomial: $(CHECK_DIR)/check_lbinomial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lbinomial $(CHECK_DIR)/check_lbinomial.cc -lquadmath

$(CHECK_DIR)/check_ldouble_factorial: $(CHECK_DIR)/check_ldouble_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_ldouble_factorial $(CHECK_DIR)/check_ldouble_factorial.cc -lquadmath

$(CHECK_DIR)/check_legendre: $(CHECK_DIR)/check_legendre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_legendre $(CHECK_DIR)/check_legendre.cc -lquadmath

$(CHECK_DIR)/check_legendre_q: $(CHECK_DIR)/check_legendre_q.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_legendre_q $(CHECK_DIR)/check_legendre_q.cc -lquadmath

$(CHECK_DIR)/check_lfactorial: $(CHECK_DIR)/check_lfactorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lfactorial $(CHECK_DIR)/check_lfactorial.cc -lquadmath

$(CHECK_DIR)/check_lfalling_factorial: $(CHECK_DIR)/check_lfalling_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lfalling_factorial $(CHECK_DIR)/check_lfalling_factorial.cc -lquadmath

$(CHECK_DIR)/check_lgamma: $(CHECK_DIR)/check_lgamma.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lgamma $(CHECK_DIR)/check_lgamma.cc -lquadmath

$(CHECK_DIR)/check_logistic_p: $(CHECK_DIR)/check_logistic_p.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_logistic_p $(CHECK_DIR)/check_logistic_p.cc -lquadmath

$(CHECK_DIR)/check_logistic_pdf: $(CHECK_DIR)/check_logistic_pdf.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_logistic_pdf $(CHECK_DIR)/check_logistic_pdf.cc -lquadmath

$(CHECK_DIR)/check_lognormal_p: $(CHECK_DIR)/check_lognormal_p.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lognormal_p $(CHECK_DIR)/check_lognormal_p.cc -lquadmath

$(CHECK_DIR)/check_lognormal_pdf: $(CHECK_DIR)/check_lognormal_pdf.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lognormal_pdf $(CHECK_DIR)/check_lognormal_pdf.cc -lquadmath

$(CHECK_DIR)/check_lrising_factorial: $(CHECK_DIR)/check_lrising_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_lrising_factorial $(CHECK_DIR)/check_lrising_factorial.cc -lquadmath

$(CHECK_DIR)/check_normal_p: $(CHECK_DIR)/check_normal_p.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_normal_p $(CHECK_DIR)/check_normal_p.cc -lquadmath

$(CHECK_DIR)/check_normal_pdf: $(CHECK_DIR)/check_normal_pdf.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_normal_pdf $(CHECK_DIR)/check_normal_pdf.cc -lquadmath

$(CHECK_DIR)/check_owens_t: $(CHECK_DIR)/check_owens_t.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_owens_t $(CHECK_DIR)/check_owens_t.cc -lquadmath

$(CHECK_DIR)/check_polygamma: $(CHECK_DIR)/check_polygamma.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_polygamma $(CHECK_DIR)/check_polygamma.cc -lquadmath

$(CHECK_DIR)/check_radpoly: $(CHECK_DIR)/check_radpoly.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_radpoly $(CHECK_DIR)/check_radpoly.cc -lquadmath

$(CHECK_DIR)/check_riemann_zeta: $(CHECK_DIR)/check_riemann_zeta.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_riemann_zeta $(CHECK_DIR)/check_riemann_zeta.cc -lquadmath

$(CHECK_DIR)/check_rising_factorial: $(CHECK_DIR)/check_rising_factorial.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_rising_factorial $(CHECK_DIR)/check_rising_factorial.cc -lquadmath

$(CHECK_DIR)/check_shi: $(CHECK_DIR)/check_shi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_shi $(CHECK_DIR)/check_shi.cc -lquadmath

$(CHECK_DIR)/check_sinc: $(CHECK_DIR)/check_sinc.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sinc $(CHECK_DIR)/check_sinc.cc -lquadmath

$(CHECK_DIR)/check_sinc_pi: $(CHECK_DIR)/check_sinc_pi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sinc_pi $(CHECK_DIR)/check_sinc_pi.cc -lquadmath

$(CHECK_DIR)/check_sinhint: $(CHECK_DIR)/check_sinhint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sinhint $(CHECK_DIR)/check_sinhint.cc -lquadmath

$(CHECK_DIR)/check_sinint: $(CHECK_DIR)/check_sinint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sinint $(CHECK_DIR)/check_sinint.cc -lquadmath

$(CHECK_DIR)/check_sin_pi: $(CHECK_DIR)/check_sin_pi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sin_pi $(CHECK_DIR)/check_sin_pi.cc -lquadmath

$(CHECK_DIR)/check_sph_bessel: $(CHECK_DIR)/check_sph_bessel.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_bessel $(CHECK_DIR)/check_sph_bessel.cc -lquadmath

$(CHECK_DIR)/check_sph_bessel_i: $(CHECK_DIR)/check_sph_bessel_i.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_bessel_i $(CHECK_DIR)/check_sph_bessel_i.cc -lquadmath

$(CHECK_DIR)/check_sph_bessel_k: $(CHECK_DIR)/check_sph_bessel_k.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_bessel_k $(CHECK_DIR)/check_sph_bessel_k.cc -lquadmath

$(CHECK_DIR)/check_sph_hankel_1: $(CHECK_DIR)/check_sph_hankel_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_hankel_1 $(CHECK_DIR)/check_sph_hankel_1.cc -lquadmath

$(CHECK_DIR)/check_sph_hankel_2: $(CHECK_DIR)/check_sph_hankel_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_hankel_2 $(CHECK_DIR)/check_sph_hankel_2.cc -lquadmath

$(CHECK_DIR)/check_sph_harmonic: $(CHECK_DIR)/check_sph_harmonic.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_harmonic $(CHECK_DIR)/check_sph_harmonic.cc -lquadmath

$(CHECK_DIR)/check_sph_legendre: $(CHECK_DIR)/check_sph_legendre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_legendre $(CHECK_DIR)/check_sph_legendre.cc -lquadmath

$(CHECK_DIR)/check_sph_neumann: $(CHECK_DIR)/check_sph_neumann.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_sph_neumann $(CHECK_DIR)/check_sph_neumann.cc -lquadmath

$(CHECK_DIR)/check_stirling_1: $(CHECK_DIR)/check_stirling_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_stirling_1 $(CHECK_DIR)/check_stirling_1.cc -lquadmath

$(CHECK_DIR)/check_stirling_2: $(CHECK_DIR)/check_stirling_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_stirling_2 $(CHECK_DIR)/check_stirling_2.cc -lquadmath

$(CHECK_DIR)/check_tgamma_lower: $(CHECK_DIR)/check_tgamma_lower.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tgamma_lower $(CHECK_DIR)/check_tgamma_lower.cc -lquadmath

$(CHECK_DIR)/check_tgamma: $(CHECK_DIR)/check_tgamma.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tgamma $(CHECK_DIR)/check_tgamma.cc -lquadmath

$(CHECK_DIR)/check_theta_1: $(CHECK_DIR)/check_theta_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_theta_1 $(CHECK_DIR)/check_theta_1.cc -lquadmath

$(CHECK_DIR)/check_theta_2: $(CHECK_DIR)/check_theta_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_theta_2 $(CHECK_DIR)/check_theta_2.cc -lquadmath

$(CHECK_DIR)/check_theta_3: $(CHECK_DIR)/check_theta_3.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_theta_3 $(CHECK_DIR)/check_theta_3.cc -lquadmath

$(CHECK_DIR)/check_theta_4: $(CHECK_DIR)/check_theta_4.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_theta_4 $(CHECK_DIR)/check_theta_4.cc -lquadmath

$(CHECK_DIR)/check_zernike: $(CHECK_DIR)/check_zernike.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_zernike $(CHECK_DIR)/check_zernike.cc -lquadmath

$(CHECK_DIR)/complex_ellint_rc: $(CHECK_DIR)/complex_ellint_rc.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_ellint_rc $(CHECK_DIR)/complex_ellint_rc.cc -lquadmath

$(CHECK_DIR)/complex_ellint_rd: $(CHECK_DIR)/complex_ellint_rd.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_ellint_rd $(CHECK_DIR)/complex_ellint_rd.cc -lquadmath

$(CHECK_DIR)/complex_ellint_rf: $(CHECK_DIR)/complex_ellint_rf.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_ellint_rf $(CHECK_DIR)/complex_ellint_rf.cc -lquadmath

$(CHECK_DIR)/complex_ellint_rg: $(CHECK_DIR)/complex_ellint_rg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_ellint_rg $(CHECK_DIR)/complex_ellint_rg.cc -lquadmath

$(CHECK_DIR)/complex_ellint_rj: $(CHECK_DIR)/complex_ellint_rj.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_ellint_rj $(CHECK_DIR)/complex_ellint_rj.cc -lquadmath

$(CHECK_DIR)/complex_airy_ai: $(CHECK_DIR)/complex_airy_ai.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_airy_ai $(CHECK_DIR)/complex_airy_ai.cc -lquadmath

$(CHECK_DIR)/complex_airy_bi: $(CHECK_DIR)/complex_airy_bi.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/complex_airy_bi $(CHECK_DIR)/complex_airy_bi.cc -lquadmath

$(CHECK_DIR)/deathmatch_comp_ellint: $(CHECK_DIR)/deathmatch_comp_ellint.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/deathmatch_comp_ellint $(CHECK_DIR)/deathmatch_comp_ellint.cc -lquadmath

$(CHECK_DIR)/deathmatch_conf_hyperg: $(CHECK_DIR)/deathmatch_conf_hyperg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/deathmatch_conf_hyperg $(CHECK_DIR)/deathmatch_conf_hyperg.cc -lquadmath

$(CHECK_DIR)/deathmatch_conf_hyperg_lim: $(CHECK_DIR)/deathmatch_conf_hyperg_lim.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/deathmatch_conf_hyperg_lim $(CHECK_DIR)/deathmatch_conf_hyperg_lim.cc -lquadmath

$(CHECK_DIR)/deathmatch_hyperg: $(CHECK_DIR)/deathmatch_hyperg.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/deathmatch_hyperg $(CHECK_DIR)/deathmatch_hyperg.cc -lquadmath

$(CHECK_DIR)/pr56216_cyl_hankel_1: $(CHECK_DIR)/pr56216_cyl_hankel_1.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr56216_cyl_hankel_1 $(CHECK_DIR)/pr56216_cyl_hankel_1.cc -lquadmath

$(CHECK_DIR)/pr56216_cyl_hankel_2: $(CHECK_DIR)/pr56216_cyl_hankel_2.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr56216_cyl_hankel_2 $(CHECK_DIR)/pr56216_cyl_hankel_2.cc -lquadmath

$(CHECK_DIR)/pr56216_cyl_bessel_i: $(CHECK_DIR)/pr56216_cyl_bessel_i.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr56216_cyl_bessel_i $(CHECK_DIR)/pr56216_cyl_bessel_i.cc -lquadmath

$(CHECK_DIR)/pr86655_assoc_legendre: $(CHECK_DIR)/pr86655_assoc_legendre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr86655_assoc_legendre $(CHECK_DIR)/pr86655_assoc_legendre.cc -lquadmath

$(CHECK_DIR)/pr86655_sph_legendre: $(CHECK_DIR)/pr86655_sph_legendre.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr86655_sph_legendre $(CHECK_DIR)/pr86655_sph_legendre.cc -lquadmath

$(CHECK_DIR)/pr68397: $(CHECK_DIR)/pr68397.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr68397 $(CHECK_DIR)/pr68397.cc -lquadmath

$(CHECK_DIR)/origin_cyl_bessel_j: $(CHECK_DIR)/origin_cyl_bessel_j.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/origin_cyl_bessel_j $(CHECK_DIR)/origin_cyl_bessel_j.cc -lquadmath

$(CHECK_DIR)/origin_cyl_neumann: $(CHECK_DIR)/origin_cyl_neumann.cc
	$(CXX17) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/origin_cyl_neumann $(CHECK_DIR)/origin_cyl_neumann.cc -lquadmath

$(CHECK_DIR)/check_tr1_assoc_laguerre: $(CHECK_DIR)/check_tr1_assoc_laguerre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_assoc_laguerre $(CHECK_DIR)/check_tr1_assoc_laguerre.cc -lquadmath

$(CHECK_DIR)/check_tr1_assoc_legendre: $(CHECK_DIR)/check_tr1_assoc_legendre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_assoc_legendre $(CHECK_DIR)/check_tr1_assoc_legendre.cc -lquadmath

$(CHECK_DIR)/check_tr1_beta: $(CHECK_DIR)/check_tr1_beta.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_beta $(CHECK_DIR)/check_tr1_beta.cc -lquadmath

$(CHECK_DIR)/check_tr1_comp_ellint_1: $(CHECK_DIR)/check_tr1_comp_ellint_1.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_comp_ellint_1 $(CHECK_DIR)/check_tr1_comp_ellint_1.cc -lquadmath

$(CHECK_DIR)/check_tr1_comp_ellint_2: $(CHECK_DIR)/check_tr1_comp_ellint_2.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_comp_ellint_2 $(CHECK_DIR)/check_tr1_comp_ellint_2.cc -lquadmath

$(CHECK_DIR)/check_tr1_comp_ellint_3: $(CHECK_DIR)/check_tr1_comp_ellint_3.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_comp_ellint_3 $(CHECK_DIR)/check_tr1_comp_ellint_3.cc -lquadmath

$(CHECK_DIR)/check_tr1_conf_hyperg: $(CHECK_DIR)/check_tr1_conf_hyperg.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_conf_hyperg $(CHECK_DIR)/check_tr1_conf_hyperg.cc -lquadmath

$(CHECK_DIR)/check_tr1_cyl_bessel_i: $(CHECK_DIR)/check_tr1_cyl_bessel_i.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_cyl_bessel_i $(CHECK_DIR)/check_tr1_cyl_bessel_i.cc -lquadmath

$(CHECK_DIR)/check_tr1_cyl_bessel_j: $(CHECK_DIR)/check_tr1_cyl_bessel_j.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_cyl_bessel_j $(CHECK_DIR)/check_tr1_cyl_bessel_j.cc -lquadmath

$(CHECK_DIR)/check_tr1_cyl_bessel_k: $(CHECK_DIR)/check_tr1_cyl_bessel_k.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_cyl_bessel_k $(CHECK_DIR)/check_tr1_cyl_bessel_k.cc -lquadmath

$(CHECK_DIR)/check_tr1_cyl_neumann: $(CHECK_DIR)/check_tr1_cyl_neumann.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_cyl_neumann $(CHECK_DIR)/check_tr1_cyl_neumann.cc -lquadmath

$(CHECK_DIR)/check_tr1_ellint_1: $(CHECK_DIR)/check_tr1_ellint_1.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_ellint_1 $(CHECK_DIR)/check_tr1_ellint_1.cc -lquadmath

$(CHECK_DIR)/check_tr1_ellint_2: $(CHECK_DIR)/check_tr1_ellint_2.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_ellint_2 $(CHECK_DIR)/check_tr1_ellint_2.cc -lquadmath

$(CHECK_DIR)/check_tr1_ellint_3: $(CHECK_DIR)/check_tr1_ellint_3.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_ellint_3 $(CHECK_DIR)/check_tr1_ellint_3.cc -lquadmath

$(CHECK_DIR)/check_tr1_expint: $(CHECK_DIR)/check_tr1_expint.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_expint $(CHECK_DIR)/check_tr1_expint.cc -lquadmath

$(CHECK_DIR)/check_tr1_hermite: $(CHECK_DIR)/check_tr1_hermite.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_hermite $(CHECK_DIR)/check_tr1_hermite.cc -lquadmath

$(CHECK_DIR)/check_tr1_hyperg: $(CHECK_DIR)/check_tr1_hyperg.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_hyperg $(CHECK_DIR)/check_tr1_hyperg.cc -lquadmath

$(CHECK_DIR)/check_tr1_laguerre: $(CHECK_DIR)/check_tr1_laguerre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_laguerre $(CHECK_DIR)/check_tr1_laguerre.cc -lquadmath

$(CHECK_DIR)/check_tr1_legendre: $(CHECK_DIR)/check_tr1_legendre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_legendre $(CHECK_DIR)/check_tr1_legendre.cc -lquadmath

$(CHECK_DIR)/check_tr1_riemann_zeta: $(CHECK_DIR)/check_tr1_riemann_zeta.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_riemann_zeta $(CHECK_DIR)/check_tr1_riemann_zeta.cc -lquadmath

$(CHECK_DIR)/check_tr1_sph_bessel: $(CHECK_DIR)/check_tr1_sph_bessel.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_sph_bessel $(CHECK_DIR)/check_tr1_sph_bessel.cc -lquadmath

$(CHECK_DIR)/check_tr1_sph_legendre: $(CHECK_DIR)/check_tr1_sph_legendre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_sph_legendre $(CHECK_DIR)/check_tr1_sph_legendre.cc -lquadmath

$(CHECK_DIR)/check_tr1_sph_neumann: $(CHECK_DIR)/check_tr1_sph_neumann.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/check_tr1_sph_neumann $(CHECK_DIR)/check_tr1_sph_neumann.cc -lquadmath

$(CHECK_DIR)/pr86655_tr1_assoc_legendre: $(CHECK_DIR)/pr86655_tr1_assoc_legendre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr86655_tr1_assoc_legendre $(CHECK_DIR)/pr86655_tr1_assoc_legendre.cc -lquadmath

$(CHECK_DIR)/pr86655_tr1_sph_legendre: $(CHECK_DIR)/pr86655_tr1_sph_legendre.cc
	$(CXX14) $(INCLUDES) -I$(CXX_TEST_INC_DIR) -D_GLIBCXX_ASSERT -D__TEST_DEBUG -o $(CHECK_DIR)/pr86655_tr1_sph_legendre $(CHECK_DIR)/pr86655_tr1_sph_legendre.cc -lquadmath


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
