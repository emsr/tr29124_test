cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_special_functions
  VERSION 0.9
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_library(cxx_special_functions INTERFACE)
target_include_directories(cxx_special_functions INTERFACE include)
target_link_libraries(cxx_special_functions INTERFACE cxx_complex_math cxx_numeric_limits)

add_library(test_support INTERFACE)
target_include_directories(test_support INTERFACE test/include)

add_executable(testcase test/testcase.cpp)
target_link_libraries(testcase wrap_gsl wrap_boost wrap_burkhardt test_support)
target_include_directories(testcase PRIVATE .)

# Function checks...

set(EMSR_LIBRARIES
  cxx_numeric_limits
  cxx_special_functions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_continued_fractions
  cxx_summation
  cxx_polynomial
  cxx_integration
  test_support)

add_executable(check_airy_ai test/src/check_airy_ai.cc)
target_include_directories(check_airy_ai PRIVATE .)
target_link_libraries(check_airy_ai ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_airy_ai COMMAND check_airy_ai)

add_executable(check_airy_ai_scaled test/src/check_airy_ai_scaled.cc)
target_include_directories(check_airy_ai_scaled PRIVATE .)
target_link_libraries(check_airy_ai_scaled ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_airy_ai_scaled COMMAND check_airy_ai_scaled)

add_executable(check_airy_bi test/src/check_airy_bi.cc)
target_include_directories(check_airy_bi PRIVATE .)
target_link_libraries(check_airy_bi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_airy_bi COMMAND check_airy_bi)

add_executable(check_airy_bi_scaled test/src/check_airy_bi_scaled.cc)
target_include_directories(check_airy_bi_scaled PRIVATE .)
target_link_libraries(check_airy_bi_scaled ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_airy_bi_scaled COMMAND check_airy_bi_scaled)

add_executable(check_assoc_laguerre test/src/check_assoc_laguerre.cc)
target_include_directories(check_assoc_laguerre PRIVATE .)
target_link_libraries(check_assoc_laguerre ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_assoc_laguerre COMMAND check_assoc_laguerre)

add_executable(check_assoc_legendre test/src/check_assoc_legendre.cc)
target_include_directories(check_assoc_legendre PRIVATE .)
target_link_libraries(check_assoc_legendre ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_assoc_legendre COMMAND check_assoc_legendre)

add_executable(check_bell test/src/check_bell.cc)
target_include_directories(check_bell PRIVATE .)
target_link_libraries(check_bell ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_bell COMMAND check_bell)

add_executable(check_bernoulli test/src/check_bernoulli.cc)
target_include_directories(check_bernoulli PRIVATE .)
target_link_libraries(check_bernoulli ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_bernoulli COMMAND check_bernoulli)

add_executable(check_beta test/src/check_beta.cc)
target_include_directories(check_beta PRIVATE .)
target_link_libraries(check_beta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_beta COMMAND check_beta)

add_executable(check_binomial test/src/check_binomial.cc)
target_include_directories(check_binomial PRIVATE .)
target_link_libraries(check_binomial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_binomial COMMAND check_binomial)

add_executable(check_chebyshev_t test/src/check_chebyshev_t.cc)
target_include_directories(check_chebyshev_t PRIVATE .)
target_link_libraries(check_chebyshev_t ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_chebyshev_t COMMAND check_chebyshev_t)

add_executable(check_chebyshev_u test/src/check_chebyshev_u.cc)
target_include_directories(check_chebyshev_u PRIVATE .)
target_link_libraries(check_chebyshev_u ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_chebyshev_u COMMAND check_chebyshev_u)

add_executable(check_chebyshev_v test/src/check_chebyshev_v.cc)
target_include_directories(check_chebyshev_v PRIVATE .)
target_link_libraries(check_chebyshev_v ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_chebyshev_v COMMAND check_chebyshev_v)

add_executable(check_chebyshev_w test/src/check_chebyshev_w.cc)
target_include_directories(check_chebyshev_w PRIVATE .)
target_link_libraries(check_chebyshev_w ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_chebyshev_w COMMAND check_chebyshev_w)

add_executable(check_chi test/src/check_chi.cc)
target_include_directories(check_chi PRIVATE .)
target_link_libraries(check_chi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_chi COMMAND check_chi)

add_executable(check_clausen_cl test/src/check_clausen_cl.cc)
target_include_directories(check_clausen_cl PRIVATE .)
target_link_libraries(check_clausen_cl ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_clausen_cl COMMAND check_clausen_cl)

add_executable(check_comp_ellint_1 test/src/check_comp_ellint_1.cc)
target_include_directories(check_comp_ellint_1 PRIVATE .)
target_link_libraries(check_comp_ellint_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_comp_ellint_1 COMMAND check_comp_ellint_1)

add_executable(check_comp_ellint_2 test/src/check_comp_ellint_2.cc)
target_include_directories(check_comp_ellint_2 PRIVATE .)
target_link_libraries(check_comp_ellint_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_comp_ellint_2 COMMAND check_comp_ellint_2)

add_executable(check_comp_ellint_3 test/src/check_comp_ellint_3.cc)
target_include_directories(check_comp_ellint_3 PRIVATE .)
target_link_libraries(check_comp_ellint_3 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_comp_ellint_3 COMMAND check_comp_ellint_3)

add_executable(check_comp_ellint_d test/src/check_comp_ellint_d.cc)
target_include_directories(check_comp_ellint_d PRIVATE .)
target_link_libraries(check_comp_ellint_d ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_comp_ellint_d COMMAND check_comp_ellint_d)

add_executable(check_conf_hyperg test/src/check_conf_hyperg.cc)
target_include_directories(check_conf_hyperg PRIVATE .)
target_link_libraries(check_conf_hyperg ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_conf_hyperg COMMAND check_conf_hyperg)

add_executable(check_conf_hyperg_lim test/src/check_conf_hyperg_lim.cc)
target_include_directories(check_conf_hyperg_lim PRIVATE .)
target_link_libraries(check_conf_hyperg_lim ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_conf_hyperg_lim COMMAND check_conf_hyperg_lim)

add_executable(check_coshint test/src/check_coshint.cc)
target_include_directories(check_coshint PRIVATE .)
target_link_libraries(check_coshint ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_coshint COMMAND check_coshint)

add_executable(check_cosint test/src/check_cosint.cc)
target_include_directories(check_cosint PRIVATE .)
target_link_libraries(check_cosint ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cosint COMMAND check_cosint)

add_executable(check_cos_pi test/src/check_cos_pi.cc)
target_include_directories(check_cos_pi PRIVATE .)
target_link_libraries(check_cos_pi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cos_pi COMMAND check_cos_pi)

add_executable(check_coulomb_f test/src/check_coulomb_f.cc)
target_include_directories(check_coulomb_f PRIVATE .)
target_link_libraries(check_coulomb_f ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_coulomb_f COMMAND check_coulomb_f)

add_executable(check_coulomb_g test/src/check_coulomb_g.cc)
target_include_directories(check_coulomb_g PRIVATE .)
target_link_libraries(check_coulomb_g ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_coulomb_g COMMAND check_coulomb_g)

add_executable(check_cyl_bessel_i test/src/check_cyl_bessel_i.cc)
target_include_directories(check_cyl_bessel_i PRIVATE .)
target_link_libraries(check_cyl_bessel_i ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_bessel_i COMMAND check_cyl_bessel_i)

add_executable(check_cyl_bessel_i_scaled test/src/check_cyl_bessel_i_scaled.cc)
target_include_directories(check_cyl_bessel_i_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_i_scaled ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_bessel_i_scaled COMMAND check_cyl_bessel_i_scaled)

add_executable(check_cyl_bessel_j test/src/check_cyl_bessel_j.cc)
target_include_directories(check_cyl_bessel_j PRIVATE .)
target_link_libraries(check_cyl_bessel_j ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_bessel_j COMMAND check_cyl_bessel_j)

add_executable(check_cyl_bessel_k test/src/check_cyl_bessel_k.cc)
target_include_directories(check_cyl_bessel_k PRIVATE .)
target_link_libraries(check_cyl_bessel_k ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_bessel_k COMMAND check_cyl_bessel_k)

add_executable(check_cyl_bessel_k_scaled test/src/check_cyl_bessel_k_scaled.cc)
target_include_directories(check_cyl_bessel_k_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_k_scaled ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_bessel_k_scaled COMMAND check_cyl_bessel_k_scaled)

add_executable(check_cyl_hankel_1 test/src/check_cyl_hankel_1.cc)
target_include_directories(check_cyl_hankel_1 PRIVATE .)
target_link_libraries(check_cyl_hankel_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_hankel_1 COMMAND check_cyl_hankel_1)

add_executable(check_cyl_hankel_2 test/src/check_cyl_hankel_2.cc)
target_include_directories(check_cyl_hankel_2 PRIVATE .)
target_link_libraries(check_cyl_hankel_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_hankel_2 COMMAND check_cyl_hankel_2)

add_executable(check_cyl_neumann test/src/check_cyl_neumann.cc)
target_include_directories(check_cyl_neumann PRIVATE .)
target_link_libraries(check_cyl_neumann ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_cyl_neumann COMMAND check_cyl_neumann)

add_executable(check_dawson test/src/check_dawson.cc)
target_include_directories(check_dawson PRIVATE .)
target_link_libraries(check_dawson ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_dawson COMMAND check_dawson)

add_executable(check_debye test/src/check_debye.cc)
target_include_directories(check_debye PRIVATE .)
target_link_libraries(check_debye ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_debye COMMAND check_debye)

add_executable(check_digamma test/src/check_digamma.cc)
target_include_directories(check_digamma PRIVATE .)
target_link_libraries(check_digamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_digamma COMMAND check_digamma)

add_executable(check_dilog test/src/check_dilog.cc)
target_include_directories(check_dilog PRIVATE .)
target_link_libraries(check_dilog ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_dilog COMMAND check_dilog)

add_executable(check_dirichlet_beta test/src/check_dirichlet_beta.cc)
target_include_directories(check_dirichlet_beta PRIVATE .)
target_link_libraries(check_dirichlet_beta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_dirichlet_beta COMMAND check_dirichlet_beta)

add_executable(check_dirichlet_eta test/src/check_dirichlet_eta.cc)
target_include_directories(check_dirichlet_eta PRIVATE .)
target_link_libraries(check_dirichlet_eta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_dirichlet_eta COMMAND check_dirichlet_eta)

add_executable(check_dirichlet_lambda test/src/check_dirichlet_lambda.cc)
target_include_directories(check_dirichlet_lambda PRIVATE .)
target_link_libraries(check_dirichlet_lambda ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_dirichlet_lambda COMMAND check_dirichlet_lambda)

add_executable(check_double_factorial test/src/check_double_factorial.cc)
target_include_directories(check_double_factorial PRIVATE .)
target_link_libraries(check_double_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_double_factorial COMMAND check_double_factorial)

add_executable(check_ellint_1 test/src/check_ellint_1.cc)
target_include_directories(check_ellint_1 PRIVATE .)
target_link_libraries(check_ellint_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_1 COMMAND check_ellint_1)

add_executable(check_ellint_2 test/src/check_ellint_2.cc)
target_include_directories(check_ellint_2 PRIVATE .)
target_link_libraries(check_ellint_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_2 COMMAND check_ellint_2)

add_executable(check_ellint_3 test/src/check_ellint_3.cc)
target_include_directories(check_ellint_3 PRIVATE .)
target_link_libraries(check_ellint_3 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_3 COMMAND check_ellint_3)

add_executable(check_ellint_d test/src/check_ellint_d.cc)
target_include_directories(check_ellint_d PRIVATE .)
target_link_libraries(check_ellint_d ${EMSR_LIBRARIES} quadmath)

add_executable(check_ellint_rc test/src/check_ellint_rc.cc)
target_include_directories(check_ellint_rc PRIVATE .)
target_link_libraries(check_ellint_rc ${EMSR_LIBRARIES} quadmath)

add_executable(check_ellint_rd test/src/check_ellint_rd.cc)
target_include_directories(check_ellint_rd PRIVATE .)
target_link_libraries(check_ellint_rd ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_rd COMMAND check_ellint_rd)

add_executable(check_ellint_rf test/src/check_ellint_rf.cc)
target_include_directories(check_ellint_rf PRIVATE .)
target_link_libraries(check_ellint_rf ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_rf COMMAND check_ellint_rf)

add_executable(check_ellint_rg test/src/check_ellint_rg.cc)
target_include_directories(check_ellint_rg PRIVATE .)
target_link_libraries(check_ellint_rg ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_rg COMMAND check_ellint_rg)

add_executable(check_ellint_rj test/src/check_ellint_rj.cc)
target_include_directories(check_ellint_rj PRIVATE .)
target_link_libraries(check_ellint_rj ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellint_rj COMMAND check_ellint_rj)

add_executable(check_ellnome test/src/check_ellnome.cc)
target_include_directories(check_ellnome PRIVATE .)
target_link_libraries(check_ellnome ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ellnome COMMAND check_ellnome)

add_executable(check_euler test/src/check_euler.cc)
target_include_directories(check_euler PRIVATE .)
target_link_libraries(check_euler ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_euler COMMAND check_euler)

add_executable(check_eulerian_1 test/src/check_eulerian_1.cc)
target_include_directories(check_eulerian_1 PRIVATE .)
target_link_libraries(check_eulerian_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_eulerian_1 COMMAND check_eulerian_1)

add_executable(check_eulerian_2 test/src/check_eulerian_2.cc)
target_include_directories(check_eulerian_2 PRIVATE .)
target_link_libraries(check_eulerian_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_eulerian_2 COMMAND check_eulerian_2)

add_executable(check_expint test/src/check_expint.cc)
target_include_directories(check_expint PRIVATE .)
target_link_libraries(check_expint ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_expint COMMAND check_expint)

add_executable(check_expint_en test/src/check_expint_en.cc)
target_include_directories(check_expint_en PRIVATE .)
target_link_libraries(check_expint_en ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_expint_en COMMAND check_expint_en)

add_executable(check_factorial test/src/check_factorial.cc)
target_include_directories(check_factorial PRIVATE .)
target_link_libraries(check_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_factorial COMMAND check_factorial)

add_executable(check_falling_factorial test/src/check_falling_factorial.cc)
target_include_directories(check_falling_factorial PRIVATE .)
target_link_libraries(check_falling_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_falling_factorial COMMAND check_falling_factorial)

add_executable(check_fresnel_c test/src/check_fresnel_c.cc)
target_include_directories(check_fresnel_c PRIVATE .)
target_link_libraries(check_fresnel_c ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_fresnel_c COMMAND check_fresnel_c)

add_executable(check_fresnel_s test/src/check_fresnel_s.cc)
target_include_directories(check_fresnel_s PRIVATE .)
target_link_libraries(check_fresnel_s ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_fresnel_s COMMAND check_fresnel_s)

add_executable(check_gamma_p test/src/check_gamma_p.cc)
target_include_directories(check_gamma_p PRIVATE .)
target_link_libraries(check_gamma_p ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_gamma_p COMMAND check_gamma_p)

add_executable(check_gamma_q test/src/check_gamma_q.cc)
target_include_directories(check_gamma_q PRIVATE .)
target_link_libraries(check_gamma_q ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_gamma_q COMMAND check_gamma_q)

add_executable(check_gamma_reciprocal test/src/check_gamma_reciprocal.cc)
target_include_directories(check_gamma_reciprocal PRIVATE .)
target_link_libraries(check_gamma_reciprocal ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_gamma_reciprocal COMMAND check_gamma_reciprocal)

add_executable(check_gegenbauer test/src/check_gegenbauer.cc)
target_include_directories(check_gegenbauer PRIVATE .)
target_link_libraries(check_gegenbauer ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_gegenbauer COMMAND check_gegenbauer)

add_executable(check_hermite test/src/check_hermite.cc)
target_include_directories(check_hermite PRIVATE .)
target_link_libraries(check_hermite ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_hermite COMMAND check_hermite)

add_executable(check_heuman_lambda test/src/check_heuman_lambda.cc)
target_include_directories(check_heuman_lambda PRIVATE .)
target_link_libraries(check_heuman_lambda ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_heuman_lambda COMMAND check_heuman_lambda)

add_executable(check_hurwitz_zeta test/src/check_hurwitz_zeta.cc)
target_include_directories(check_hurwitz_zeta PRIVATE .)
target_link_libraries(check_hurwitz_zeta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_hurwitz_zeta COMMAND check_hurwitz_zeta)

add_executable(check_hyperg test/src/check_hyperg.cc)
target_include_directories(check_hyperg PRIVATE .)
target_link_libraries(check_hyperg ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_hyperg COMMAND check_hyperg)

add_executable(check_ibeta test/src/check_ibeta.cc)
target_include_directories(check_ibeta PRIVATE .)
target_link_libraries(check_ibeta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ibeta COMMAND check_ibeta)

add_executable(check_ibetac test/src/check_ibetac.cc)
target_include_directories(check_ibetac PRIVATE .)
target_link_libraries(check_ibetac ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ibetac COMMAND check_ibetac)

add_executable(check_jacobi test/src/check_jacobi.cc)
target_include_directories(check_jacobi PRIVATE .)
target_link_libraries(check_jacobi ${EMSR_LIBRARIES} quadmath)

add_executable(check_jacobi_cn test/src/check_jacobi_cn.cc)
target_include_directories(check_jacobi_cn PRIVATE .)
target_link_libraries(check_jacobi_cn ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_jacobi_cn COMMAND check_jacobi_cn)

add_executable(check_jacobi_dn test/src/check_jacobi_dn.cc)
target_include_directories(check_jacobi_dn PRIVATE .)
target_link_libraries(check_jacobi_dn ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_jacobi_dn COMMAND check_jacobi_dn)

add_executable(check_jacobi_sn test/src/check_jacobi_sn.cc)
target_include_directories(check_jacobi_sn PRIVATE .)
target_link_libraries(check_jacobi_sn ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_jacobi_sn COMMAND check_jacobi_sn)

add_executable(check_jacobi_zeta test/src/check_jacobi_zeta.cc)
target_include_directories(check_jacobi_zeta PRIVATE .)
target_link_libraries(check_jacobi_zeta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_jacobi_zeta COMMAND check_jacobi_zeta)

add_executable(check_laguerre test/src/check_laguerre.cc)
target_include_directories(check_laguerre PRIVATE .)
target_link_libraries(check_laguerre ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_laguerre COMMAND check_laguerre)

add_executable(check_lah test/src/check_lah.cc)
target_include_directories(check_lah PRIVATE .)
target_link_libraries(check_lah ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lah COMMAND check_lah)

add_executable(check_lbinomial test/src/check_lbinomial.cc)
target_include_directories(check_lbinomial PRIVATE .)
target_link_libraries(check_lbinomial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lbinomial COMMAND check_lbinomial)

add_executable(check_ldouble_factorial test/src/check_ldouble_factorial.cc)
target_include_directories(check_ldouble_factorial PRIVATE .)
target_link_libraries(check_ldouble_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_ldouble_factorial COMMAND check_ldouble_factorial)

add_executable(check_legendre test/src/check_legendre.cc)
target_include_directories(check_legendre PRIVATE .)
target_link_libraries(check_legendre ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_legendre COMMAND check_legendre)

add_executable(check_legendre_q test/src/check_legendre_q.cc)
target_include_directories(check_legendre_q PRIVATE .)
target_link_libraries(check_legendre_q ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_legendre_q COMMAND check_legendre_q)

add_executable(check_lfactorial test/src/check_lfactorial.cc)
target_include_directories(check_lfactorial PRIVATE .)
target_link_libraries(check_lfactorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lfactorial COMMAND check_lfactorial)

add_executable(check_lfalling_factorial test/src/check_lfalling_factorial.cc)
target_include_directories(check_lfalling_factorial PRIVATE .)
target_link_libraries(check_lfalling_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lfalling_factorial COMMAND check_lfalling_factorial)

add_executable(check_lgamma test/src/check_lgamma.cc)
target_include_directories(check_lgamma PRIVATE .)
target_link_libraries(check_lgamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lgamma COMMAND check_lgamma)

add_executable(check_logistic_p test/src/check_logistic_p.cc)
target_include_directories(check_logistic_p PRIVATE .)
target_link_libraries(check_logistic_p ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_logistic_p COMMAND check_logistic_p)

add_executable(check_logistic_pdf test/src/check_logistic_pdf.cc)
target_include_directories(check_logistic_pdf PRIVATE .)
target_link_libraries(check_logistic_pdf ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_logistic_pdf COMMAND check_logistic_pdf)

add_executable(check_lognormal_p test/src/check_lognormal_p.cc)
target_include_directories(check_lognormal_p PRIVATE .)
target_link_libraries(check_lognormal_p ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lognormal_p COMMAND check_lognormal_p)

add_executable(check_lognormal_pdf test/src/check_lognormal_pdf.cc)
target_include_directories(check_lognormal_pdf PRIVATE .)
target_link_libraries(check_lognormal_pdf ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lognormal_pdf COMMAND check_lognormal_pdf)

add_executable(check_lrising_factorial test/src/check_lrising_factorial.cc)
target_include_directories(check_lrising_factorial PRIVATE .)
target_link_libraries(check_lrising_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_lrising_factorial COMMAND check_lrising_factorial)

add_executable(check_normal_p test/src/check_normal_p.cc)
target_include_directories(check_normal_p PRIVATE .)
target_link_libraries(check_normal_p ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_normal_p COMMAND check_normal_p)

add_executable(check_normal_pdf test/src/check_normal_pdf.cc)
target_include_directories(check_normal_pdf PRIVATE .)
target_link_libraries(check_normal_pdf ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_normal_pdf COMMAND check_normal_pdf)

add_executable(check_owens_t test/src/check_owens_t.cc)
target_include_directories(check_owens_t PRIVATE .)
target_link_libraries(check_owens_t ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_owens_t COMMAND check_owens_t)

add_executable(check_polygamma test/src/check_polygamma.cc)
target_include_directories(check_polygamma PRIVATE .)
target_link_libraries(check_polygamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_polygamma COMMAND check_polygamma)

add_executable(check_radpoly test/src/check_radpoly.cc)
target_include_directories(check_radpoly PRIVATE .)
target_link_libraries(check_radpoly ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_radpoly COMMAND check_radpoly)

add_executable(check_riemann_zeta test/src/check_riemann_zeta.cc)
target_include_directories(check_riemann_zeta PRIVATE .)
target_link_libraries(check_riemann_zeta ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_riemann_zeta COMMAND check_riemann_zeta)

add_executable(check_rising_factorial test/src/check_rising_factorial.cc)
target_include_directories(check_rising_factorial PRIVATE .)
target_link_libraries(check_rising_factorial ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_rising_factorial COMMAND check_rising_factorial)

add_executable(check_shi test/src/check_shi.cc)
target_include_directories(check_shi PRIVATE .)
target_link_libraries(check_shi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_shi COMMAND check_shi)

add_executable(check_sinc test/src/check_sinc.cc)
target_include_directories(check_sinc PRIVATE .)
target_link_libraries(check_sinc ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sinc COMMAND check_sinc)

add_executable(check_sinc_pi test/src/check_sinc_pi.cc)
target_include_directories(check_sinc_pi PRIVATE .)
target_link_libraries(check_sinc_pi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sinc_pi COMMAND check_sinc_pi)

add_executable(check_sinhint test/src/check_sinhint.cc)
target_include_directories(check_sinhint PRIVATE .)
target_link_libraries(check_sinhint ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sinhint COMMAND check_sinhint)

add_executable(check_sinint test/src/check_sinint.cc)
target_include_directories(check_sinint PRIVATE .)
target_link_libraries(check_sinint ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sinint COMMAND check_sinint)

add_executable(check_sin_pi test/src/check_sin_pi.cc)
target_include_directories(check_sin_pi PRIVATE .)
target_link_libraries(check_sin_pi ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sin_pi COMMAND check_sin_pi)

add_executable(check_sph_bessel test/src/check_sph_bessel.cc)
target_include_directories(check_sph_bessel PRIVATE .)
target_link_libraries(check_sph_bessel ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_bessel COMMAND check_sph_bessel)

add_executable(check_sph_bessel_i test/src/check_sph_bessel_i.cc)
target_include_directories(check_sph_bessel_i PRIVATE .)
target_link_libraries(check_sph_bessel_i ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_bessel_i COMMAND check_sph_bessel_i)

add_executable(check_sph_bessel_k test/src/check_sph_bessel_k.cc)
target_include_directories(check_sph_bessel_k PRIVATE .)
target_link_libraries(check_sph_bessel_k ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_bessel_k COMMAND check_sph_bessel_k)

add_executable(check_sph_hankel_1 test/src/check_sph_hankel_1.cc)
target_include_directories(check_sph_hankel_1 PRIVATE .)
target_link_libraries(check_sph_hankel_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_hankel_1 COMMAND check_sph_hankel_1)

add_executable(check_sph_hankel_2 test/src/check_sph_hankel_2.cc)
target_include_directories(check_sph_hankel_2 PRIVATE .)
target_link_libraries(check_sph_hankel_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_hankel_2 COMMAND check_sph_hankel_2)

add_executable(check_sph_harmonic test/src/check_sph_harmonic.cc)
target_include_directories(check_sph_harmonic PRIVATE .)
target_link_libraries(check_sph_harmonic ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_harmonic COMMAND check_sph_harmonic)

add_executable(check_sph_legendre test/src/check_sph_legendre.cc)
target_include_directories(check_sph_legendre PRIVATE .)
target_link_libraries(check_sph_legendre ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_legendre COMMAND check_sph_legendre)

add_executable(check_sph_neumann test/src/check_sph_neumann.cc)
target_include_directories(check_sph_neumann PRIVATE .)
target_link_libraries(check_sph_neumann ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_sph_neumann COMMAND check_sph_neumann)

add_executable(check_stirling_1 test/src/check_stirling_1.cc)
target_include_directories(check_stirling_1 PRIVATE .)
target_link_libraries(check_stirling_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_stirling_1 COMMAND check_stirling_1)

add_executable(check_stirling_2 test/src/check_stirling_2.cc)
target_include_directories(check_stirling_2 PRIVATE .)
target_link_libraries(check_stirling_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_stirling_2 COMMAND check_stirling_2)

add_executable(check_tgamma test/src/check_tgamma.cc)
target_include_directories(check_tgamma PRIVATE .)
target_link_libraries(check_tgamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_tgamma COMMAND check_tgamma)

add_executable(check_tgamma_lower test/src/check_tgamma_lower.cc)
target_include_directories(check_tgamma_lower PRIVATE .)
target_link_libraries(check_tgamma_lower ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_tgamma_lower COMMAND check_tgamma_lower)

add_executable(check_theta_1 test/src/check_theta_1.cc)
target_include_directories(check_theta_1 PRIVATE .)
target_link_libraries(check_theta_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_theta_1 COMMAND check_theta_1)

add_executable(check_theta_2 test/src/check_theta_2.cc)
target_include_directories(check_theta_2 PRIVATE .)
target_link_libraries(check_theta_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_theta_2 COMMAND check_theta_2)

add_executable(check_theta_3 test/src/check_theta_3.cc)
target_include_directories(check_theta_3 PRIVATE .)
target_link_libraries(check_theta_3 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_theta_3 COMMAND check_theta_3)

add_executable(check_theta_4 test/src/check_theta_4.cc)
target_include_directories(check_theta_4 PRIVATE .)
target_link_libraries(check_theta_4 ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_theta_4 COMMAND check_theta_4)

add_executable(check_zernike test/src/check_zernike.cc)
target_include_directories(check_zernike PRIVATE .)
target_link_libraries(check_zernike ${EMSR_LIBRARIES} quadmath)
add_test(NAME check_zernike COMMAND check_zernike)

add_executable(complex_airy_ai test/src/complex_airy_ai.cc)
target_include_directories(complex_airy_ai PRIVATE .)
target_link_libraries(complex_airy_ai ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_airy_ai COMMAND complex_airy_ai)

add_executable(complex_airy_bi test/src/complex_airy_bi.cc)
target_include_directories(complex_airy_bi PRIVATE .)
target_link_libraries(complex_airy_bi ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_airy_bi COMMAND complex_airy_bi)

add_executable(complex_ellint_rc test/src/complex_ellint_rc.cc)
target_include_directories(complex_ellint_rc PRIVATE .)
target_link_libraries(complex_ellint_rc ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_ellint_rc COMMAND complex_ellint_rc)

add_executable(complex_ellint_rd test/src/complex_ellint_rd.cc)
target_include_directories(complex_ellint_rd PRIVATE .)
target_link_libraries(complex_ellint_rd ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_ellint_rd COMMAND complex_ellint_rd)

add_executable(complex_ellint_rf test/src/complex_ellint_rf.cc)
target_include_directories(complex_ellint_rf PRIVATE .)
target_link_libraries(complex_ellint_rf ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_ellint_rf COMMAND complex_ellint_rf)

add_executable(complex_ellint_rg test/src/complex_ellint_rg.cc)
target_include_directories(complex_ellint_rg PRIVATE .)
target_link_libraries(complex_ellint_rg ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_ellint_rg COMMAND complex_ellint_rg)

add_executable(complex_ellint_rj test/src/complex_ellint_rj.cc)
target_include_directories(complex_ellint_rj PRIVATE .)
target_link_libraries(complex_ellint_rj ${EMSR_LIBRARIES} quadmath)
add_test(NAME complex_ellint_rj COMMAND complex_ellint_rj)

add_executable(deathmatch_comp_ellint test/src/deathmatch_comp_ellint.cc)
target_include_directories(deathmatch_comp_ellint PRIVATE .)
target_link_libraries(deathmatch_comp_ellint ${EMSR_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg test/src/deathmatch_conf_hyperg.cc)
target_include_directories(deathmatch_conf_hyperg PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg ${EMSR_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg_lim test/src/deathmatch_conf_hyperg_lim.cc)
target_include_directories(deathmatch_conf_hyperg_lim PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg_lim ${EMSR_LIBRARIES} quadmath)
add_test(NAME deathmatch_conf_hyperg_lim COMMAND deathmatch_conf_hyperg_lim)

add_executable(deathmatch_hyperg test/src/deathmatch_hyperg.cc)
target_include_directories(deathmatch_hyperg PRIVATE .)
target_link_libraries(deathmatch_hyperg ${EMSR_LIBRARIES} quadmath)
add_test(NAME deathmatch_hyperg COMMAND deathmatch_hyperg)

add_executable(origin_cyl_bessel_j test/src/origin_cyl_bessel_j.cc)
target_include_directories(origin_cyl_bessel_j PRIVATE .)
target_link_libraries(origin_cyl_bessel_j ${EMSR_LIBRARIES} quadmath)
add_test(NAME origin_cyl_bessel_j COMMAND origin_cyl_bessel_j)

add_executable(origin_cyl_neumann test/src/origin_cyl_neumann.cc)
target_include_directories(origin_cyl_neumann PRIVATE .)
target_link_libraries(origin_cyl_neumann ${EMSR_LIBRARIES} quadmath)
add_test(NAME origin_cyl_neumann COMMAND origin_cyl_neumann)

add_executable(pr56216_cyl_bessel_i test/src/pr56216_cyl_bessel_i.cc)
target_include_directories(pr56216_cyl_bessel_i PRIVATE .)
target_link_libraries(pr56216_cyl_bessel_i ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr56216_cyl_bessel_i COMMAND pr56216_cyl_bessel_i)

add_executable(pr56216_cyl_hankel_1 test/src/pr56216_cyl_hankel_1.cc)
target_include_directories(pr56216_cyl_hankel_1 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_1 ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr56216_cyl_hankel_1 COMMAND pr56216_cyl_hankel_1)

add_executable(pr56216_cyl_hankel_2 test/src/pr56216_cyl_hankel_2.cc)
target_include_directories(pr56216_cyl_hankel_2 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_2 ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr56216_cyl_hankel_2 COMMAND pr56216_cyl_hankel_2)

add_executable(pr68397 test/src/pr68397.cc)
target_include_directories(pr68397 PRIVATE .)
target_link_libraries(pr68397 ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr68397 COMMAND pr68397)

add_executable(pr83140 test/src/pr83140.cc)
target_include_directories(pr83140 PRIVATE .)
target_link_libraries(pr83140 ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr83140 COMMAND pr83140)

add_executable(pr86655_assoc_legendre test/src/pr86655_assoc_legendre.cc)
target_include_directories(pr86655_assoc_legendre PRIVATE .)
target_link_libraries(pr86655_assoc_legendre ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr86655_assoc_legendre COMMAND pr86655_assoc_legendre)

add_executable(pr86655_sph_legendre test/src/pr86655_sph_legendre.cc)
target_include_directories(pr86655_sph_legendre PRIVATE .)
target_link_libraries(pr86655_sph_legendre ${EMSR_LIBRARIES} quadmath)
add_test(NAME pr86655_sph_legendre COMMAND pr86655_sph_legendre)

# Orthonormality checks

add_executable(orthonorm_assoc_laguerre test/src/orthonorm_assoc_laguerre.cc)
target_include_directories(orthonorm_assoc_laguerre PRIVATE .)
target_link_libraries(orthonorm_assoc_laguerre ${EMSR_LIBRARIES})
add_test(NAME orthonorm_assoc_laguerre COMMAND orthonorm_assoc_laguerre)

add_executable(orthonorm_assoc_legendre test/src/orthonorm_assoc_legendre.cc)
target_include_directories(orthonorm_assoc_legendre PRIVATE .)
target_link_libraries(orthonorm_assoc_legendre ${EMSR_LIBRARIES})
add_test(NAME orthonorm_assoc_legendre COMMAND orthonorm_assoc_legendre)

add_executable(orthonorm_chebyshev_t test/src/orthonorm_chebyshev_t.cc)
target_include_directories(orthonorm_chebyshev_t PRIVATE .)
target_link_libraries(orthonorm_chebyshev_t ${EMSR_LIBRARIES})
add_test(NAME orthonorm_chebyshev_t COMMAND orthonorm_chebyshev_t)

add_executable(orthonorm_chebyshev_u test/src/orthonorm_chebyshev_u.cc)
target_include_directories(orthonorm_chebyshev_u PRIVATE .)
target_link_libraries(orthonorm_chebyshev_u ${EMSR_LIBRARIES})
add_test(NAME orthonorm_chebyshev_u COMMAND orthonorm_chebyshev_u)

add_executable(orthonorm_chebyshev_v test/src/orthonorm_chebyshev_v.cc)
target_include_directories(orthonorm_chebyshev_v PRIVATE .)
target_link_libraries(orthonorm_chebyshev_v ${EMSR_LIBRARIES})
add_test(NAME orthonorm_chebyshev_v COMMAND orthonorm_chebyshev_v)

add_executable(orthonorm_chebyshev_w test/src/orthonorm_chebyshev_w.cc)
target_include_directories(orthonorm_chebyshev_w PRIVATE .)
target_link_libraries(orthonorm_chebyshev_w ${EMSR_LIBRARIES})
add_test(NAME orthonorm_chebyshev_w COMMAND orthonorm_chebyshev_w)

add_executable(orthonorm_gegenbauer test/src/orthonorm_gegenbauer.cc)
target_include_directories(orthonorm_gegenbauer PRIVATE .)
target_link_libraries(orthonorm_gegenbauer ${EMSR_LIBRARIES})
add_test(NAME orthonorm_gegenbauer COMMAND orthonorm_gegenbauer)

add_executable(orthonorm_hermite test/src/orthonorm_hermite.cc)
target_include_directories(orthonorm_hermite PRIVATE .)
target_link_libraries(orthonorm_hermite ${EMSR_LIBRARIES})
add_test(NAME orthonorm_hermite COMMAND orthonorm_hermite)

add_executable(orthonorm_jacobi test/src/orthonorm_jacobi.cc)
target_include_directories(orthonorm_jacobi PRIVATE .)
target_link_libraries(orthonorm_jacobi ${EMSR_LIBRARIES})
add_test(NAME orthonorm_jacobi COMMAND orthonorm_jacobi)

add_executable(orthonorm_laguerre test/src/orthonorm_laguerre.cc)
target_include_directories(orthonorm_laguerre PRIVATE .)
target_link_libraries(orthonorm_laguerre ${EMSR_LIBRARIES})
add_test(NAME orthonorm_laguerre COMMAND orthonorm_laguerre)

add_executable(orthonorm_legendre test/src/orthonorm_legendre.cc)
target_include_directories(orthonorm_legendre PRIVATE .)
target_link_libraries(orthonorm_legendre ${EMSR_LIBRARIES})
add_test(NAME orthonorm_legendre COMMAND orthonorm_legendre)

add_executable(orthonorm_radpoly test/src/orthonorm_radpoly.cc)
target_include_directories(orthonorm_radpoly PRIVATE .)
target_link_libraries(orthonorm_radpoly ${EMSR_LIBRARIES})
add_test(NAME orthonorm_radpoly COMMAND orthonorm_radpoly)

add_executable(orthonorm_sph_legendre test/src/orthonorm_sph_legendre.cc)
target_include_directories(orthonorm_sph_legendre PRIVATE .)
target_link_libraries(orthonorm_sph_legendre ${EMSR_LIBRARIES})
add_test(NAME orthonorm_sph_legendre COMMAND orthonorm_sph_legendre)

add_executable(orthonorm_zernike test/src/orthonorm_zernike.cc)
target_include_directories(orthonorm_zernike PRIVATE .)
target_link_libraries(orthonorm_zernike ${EMSR_LIBRARIES})
add_test(NAME orthonorm_zernike COMMAND orthonorm_zernike)
