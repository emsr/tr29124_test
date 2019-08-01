cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_chebyshev
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

#add_library( INTERFACE)
#target_include_directories( INTERFACE qwerty)

add_definitions(-D__TEST_DEBUG)

set(TR29124_LIBRARIES
  cxx_special_functions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_continued_fractions
  cxx_summation
  cxx_polynomial
  cxx_integration)

add_executable(testcase testcase.cpp)
target_link_libraries(testcase wrap_gsl wrap_boost wrap_burkhardt)
target_include_directories(testcase PRIVATE .)

add_executable(testcase_tr1 testcase.cpp)
target_link_libraries(testcase_tr1 quadmath wrap_gsl wrap_boost)
target_compile_definitions(testcase_tr1 PRIVATE TR1)
target_include_directories(testcase_tr1 PRIVATE .)

#
# Function to run checks.
#
function(run_check target)
  add_test(NAME run_${target}
#    COMMAND bash -c
#      "if [ ! -f output/check_err.txt ]; then"
#      "  mkdir -p output && echo 'Check Errors' > output/check_err.txt;"
#      "fi;"
#      "if [ ! -f output/check_out.txt ]; then"
#      "  mkdir -p output && echo 'Check Output' > output/check_out.txt;"
#      "fi;"
#      "echo ${target} 2>> output/check_err.txt;"
    COMMAND bash -c
      "${CMAKE_BINARY_DIR}/bin/${target}")
endfunction()

add_executable(check_airy_ai check/check_airy_ai.cc)
target_include_directories(check_airy_ai PRIVATE .)
target_link_libraries(check_airy_ai ${TR29124_LIBRARIES} quadmath)
run_check(check_airy_ai)

add_executable(check_airy_bi check/check_airy_bi.cc)
target_include_directories(check_airy_bi PRIVATE .)
target_link_libraries(check_airy_bi ${TR29124_LIBRARIES} quadmath)
run_check(check_airy_bi)

add_executable(check_assoc_laguerre check/check_assoc_laguerre.cc)
target_include_directories(check_assoc_laguerre PRIVATE .)
target_link_libraries(check_assoc_laguerre ${TR29124_LIBRARIES} quadmath)
run_check(check_assoc_laguerre)

add_executable(check_assoc_legendre check/check_assoc_legendre.cc)
target_include_directories(check_assoc_legendre PRIVATE .)
target_link_libraries(check_assoc_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_assoc_legendre)

add_executable(check_bell check/check_bell.cc)
target_include_directories(check_bell PRIVATE .)
target_link_libraries(check_bell ${TR29124_LIBRARIES} quadmath)
run_check(check_bell)

add_executable(check_bernoulli check/check_bernoulli.cc)
target_include_directories(check_bernoulli PRIVATE .)
target_link_libraries(check_bernoulli ${TR29124_LIBRARIES} quadmath)
run_check(check_bernoulli)

add_executable(check_beta check/check_beta.cc)
target_include_directories(check_beta PRIVATE .)
target_link_libraries(check_beta ${TR29124_LIBRARIES} quadmath)
run_check(check_beta)

add_executable(check_binomial check/check_binomial.cc)
target_include_directories(check_binomial PRIVATE .)
target_link_libraries(check_binomial ${TR29124_LIBRARIES} quadmath)
run_check(check_binomial)

add_executable(check_chebyshev_t check/check_chebyshev_t.cc)
target_include_directories(check_chebyshev_t PRIVATE .)
target_link_libraries(check_chebyshev_t ${TR29124_LIBRARIES} quadmath)
run_check(check_chebyshev_t)

add_executable(check_chebyshev_u check/check_chebyshev_u.cc)
target_include_directories(check_chebyshev_u PRIVATE .)
target_link_libraries(check_chebyshev_u ${TR29124_LIBRARIES} quadmath)
run_check(check_chebyshev_u)

add_executable(check_chebyshev_v check/check_chebyshev_v.cc)
target_include_directories(check_chebyshev_v PRIVATE .)
target_link_libraries(check_chebyshev_v ${TR29124_LIBRARIES} quadmath)
run_check(check_chebyshev_v)

add_executable(check_chebyshev_w check/check_chebyshev_w.cc)
target_include_directories(check_chebyshev_w PRIVATE .)
target_link_libraries(check_chebyshev_w ${TR29124_LIBRARIES} quadmath)
run_check(check_chebyshev_w)

add_executable(check_chi check/check_chi.cc)
target_include_directories(check_chi PRIVATE .)
target_link_libraries(check_chi ${TR29124_LIBRARIES} quadmath)
run_check(check_chi)

add_executable(check_clausen_cl check/check_clausen_cl.cc)
target_include_directories(check_clausen_cl PRIVATE .)
target_link_libraries(check_clausen_cl ${TR29124_LIBRARIES} quadmath)
run_check(check_clausen_cl)

add_executable(check_comp_ellint_1 check/check_comp_ellint_1.cc)
target_include_directories(check_comp_ellint_1 PRIVATE .)
target_link_libraries(check_comp_ellint_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_comp_ellint_1)

add_executable(check_comp_ellint_2 check/check_comp_ellint_2.cc)
target_include_directories(check_comp_ellint_2 PRIVATE .)
target_link_libraries(check_comp_ellint_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_comp_ellint_2)

add_executable(check_comp_ellint_3 check/check_comp_ellint_3.cc)
target_include_directories(check_comp_ellint_3 PRIVATE .)
target_link_libraries(check_comp_ellint_3 ${TR29124_LIBRARIES} quadmath)
run_check(check_comp_ellint_3)

add_executable(check_comp_ellint_d check/check_comp_ellint_d.cc)
target_include_directories(check_comp_ellint_d PRIVATE .)
target_link_libraries(check_comp_ellint_d ${TR29124_LIBRARIES} quadmath)
run_check(check_comp_ellint_d)

add_executable(check_conf_hyperg check/check_conf_hyperg.cc)
target_include_directories(check_conf_hyperg PRIVATE .)
target_link_libraries(check_conf_hyperg ${TR29124_LIBRARIES} quadmath)
run_check(check_conf_hyperg)

add_executable(check_conf_hyperg_lim check/check_conf_hyperg_lim.cc)
target_include_directories(check_conf_hyperg_lim PRIVATE .)
target_link_libraries(check_conf_hyperg_lim ${TR29124_LIBRARIES} quadmath)
run_check(check_conf_hyperg_lim)

add_executable(check_coshint check/check_coshint.cc)
target_include_directories(check_coshint PRIVATE .)
target_link_libraries(check_coshint ${TR29124_LIBRARIES} quadmath)
run_check(check_coshint)

add_executable(check_cosint check/check_cosint.cc)
target_include_directories(check_cosint PRIVATE .)
target_link_libraries(check_cosint ${TR29124_LIBRARIES} quadmath)
run_check(check_cosint)

add_executable(check_cos_pi check/check_cos_pi.cc)
target_include_directories(check_cos_pi PRIVATE .)
target_link_libraries(check_cos_pi ${TR29124_LIBRARIES} quadmath)
run_check(check_cos_pi)

add_executable(check_cyl_bessel_i check/check_cyl_bessel_i.cc)
target_include_directories(check_cyl_bessel_i PRIVATE .)
target_link_libraries(check_cyl_bessel_i ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_bessel_i)

add_executable(check_cyl_bessel_i_scaled check/check_cyl_bessel_i_scaled.cc)
target_include_directories(check_cyl_bessel_i_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_i_scaled ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_bessel_i_scaled)

add_executable(check_cyl_bessel_j check/check_cyl_bessel_j.cc)
target_include_directories(check_cyl_bessel_j PRIVATE .)
target_link_libraries(check_cyl_bessel_j ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_bessel_j)

add_executable(check_cyl_bessel_k check/check_cyl_bessel_k.cc)
target_include_directories(check_cyl_bessel_k PRIVATE .)
target_link_libraries(check_cyl_bessel_k ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_bessel_k)

add_executable(check_cyl_bessel_k_scaled check/check_cyl_bessel_k_scaled.cc)
target_include_directories(check_cyl_bessel_k_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_k_scaled ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_bessel_k_scaled)

add_executable(check_cyl_hankel_1 check/check_cyl_hankel_1.cc)
target_include_directories(check_cyl_hankel_1 PRIVATE .)
target_link_libraries(check_cyl_hankel_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_hankel_1)

add_executable(check_cyl_hankel_2 check/check_cyl_hankel_2.cc)
target_include_directories(check_cyl_hankel_2 PRIVATE .)
target_link_libraries(check_cyl_hankel_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_hankel_2)

add_executable(check_cyl_neumann check/check_cyl_neumann.cc)
target_include_directories(check_cyl_neumann PRIVATE .)
target_link_libraries(check_cyl_neumann ${TR29124_LIBRARIES} quadmath)
run_check(check_cyl_neumann)

add_executable(check_dawson check/check_dawson.cc)
target_include_directories(check_dawson PRIVATE .)
target_link_libraries(check_dawson ${TR29124_LIBRARIES} quadmath)
run_check(check_dawson)

add_executable(check_debye check/check_debye.cc)
target_include_directories(check_debye PRIVATE .)
target_link_libraries(check_debye ${TR29124_LIBRARIES} quadmath)
run_check(check_debye)

add_executable(check_digamma check/check_digamma.cc)
target_include_directories(check_digamma PRIVATE .)
target_link_libraries(check_digamma ${TR29124_LIBRARIES} quadmath)
run_check(check_digamma)

add_executable(check_dilog check/check_dilog.cc)
target_include_directories(check_dilog PRIVATE .)
target_link_libraries(check_dilog ${TR29124_LIBRARIES} quadmath)
run_check(check_dilog)

add_executable(check_dirichlet_beta check/check_dirichlet_beta.cc)
target_include_directories(check_dirichlet_beta PRIVATE .)
target_link_libraries(check_dirichlet_beta ${TR29124_LIBRARIES} quadmath)
run_check(check_dirichlet_beta)

add_executable(check_dirichlet_eta check/check_dirichlet_eta.cc)
target_include_directories(check_dirichlet_eta PRIVATE .)
target_link_libraries(check_dirichlet_eta ${TR29124_LIBRARIES} quadmath)
run_check(check_dirichlet_eta)

add_executable(check_dirichlet_lambda check/check_dirichlet_lambda.cc)
target_include_directories(check_dirichlet_lambda PRIVATE .)
target_link_libraries(check_dirichlet_lambda ${TR29124_LIBRARIES} quadmath)
run_check(check_dirichlet_lambda)

add_executable(check_double_factorial check/check_double_factorial.cc)
target_include_directories(check_double_factorial PRIVATE .)
target_link_libraries(check_double_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_double_factorial)

add_executable(check_ellint_1 check/check_ellint_1.cc)
target_include_directories(check_ellint_1 PRIVATE .)
target_link_libraries(check_ellint_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_1)

add_executable(check_ellint_2 check/check_ellint_2.cc)
target_include_directories(check_ellint_2 PRIVATE .)
target_link_libraries(check_ellint_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_2)

add_executable(check_ellint_3 check/check_ellint_3.cc)
target_include_directories(check_ellint_3 PRIVATE .)
target_link_libraries(check_ellint_3 ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_3)

add_executable(check_ellint_d check/check_ellint_d.cc)
target_include_directories(check_ellint_d PRIVATE .)
target_link_libraries(check_ellint_d ${TR29124_LIBRARIES} quadmath)

add_executable(check_ellint_rc check/check_ellint_rc.cc)
target_include_directories(check_ellint_rc PRIVATE .)
target_link_libraries(check_ellint_rc ${TR29124_LIBRARIES} quadmath)

add_executable(check_ellint_rd check/check_ellint_rd.cc)
target_include_directories(check_ellint_rd PRIVATE .)
target_link_libraries(check_ellint_rd ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_rd)

add_executable(check_ellint_rf check/check_ellint_rf.cc)
target_include_directories(check_ellint_rf PRIVATE .)
target_link_libraries(check_ellint_rf ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_rf)

add_executable(check_ellint_rg check/check_ellint_rg.cc)
target_include_directories(check_ellint_rg PRIVATE .)
target_link_libraries(check_ellint_rg ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_rg)

add_executable(check_ellint_rj check/check_ellint_rj.cc)
target_include_directories(check_ellint_rj PRIVATE .)
target_link_libraries(check_ellint_rj ${TR29124_LIBRARIES} quadmath)
run_check(check_ellint_rj)

add_executable(check_ellnome check/check_ellnome.cc)
target_include_directories(check_ellnome PRIVATE .)
target_link_libraries(check_ellnome ${TR29124_LIBRARIES} quadmath)
run_check(check_ellnome)

add_executable(check_euler check/check_euler.cc)
target_include_directories(check_euler PRIVATE .)
target_link_libraries(check_euler ${TR29124_LIBRARIES} quadmath)
run_check(check_euler)

add_executable(check_eulerian_1 check/check_eulerian_1.cc)
target_include_directories(check_eulerian_1 PRIVATE .)
target_link_libraries(check_eulerian_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_eulerian_1)

add_executable(check_eulerian_2 check/check_eulerian_2.cc)
target_include_directories(check_eulerian_2 PRIVATE .)
target_link_libraries(check_eulerian_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_eulerian_2)

add_executable(check_expint check/check_expint.cc)
target_include_directories(check_expint PRIVATE .)
target_link_libraries(check_expint ${TR29124_LIBRARIES} quadmath)
run_check(check_expint)

add_executable(check_expint_en check/check_expint_en.cc)
target_include_directories(check_expint_en PRIVATE .)
target_link_libraries(check_expint_en ${TR29124_LIBRARIES} quadmath)
run_check(check_expint_en)

add_executable(check_factorial check/check_factorial.cc)
target_include_directories(check_factorial PRIVATE .)
target_link_libraries(check_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_factorial)

add_executable(check_falling_factorial check/check_falling_factorial.cc)
target_include_directories(check_falling_factorial PRIVATE .)
target_link_libraries(check_falling_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_falling_factorial)

add_executable(check_fresnel_c check/check_fresnel_c.cc)
target_include_directories(check_fresnel_c PRIVATE .)
target_link_libraries(check_fresnel_c ${TR29124_LIBRARIES} quadmath)
run_check(check_fresnel_c)

add_executable(check_fresnel_s check/check_fresnel_s.cc)
target_include_directories(check_fresnel_s PRIVATE .)
target_link_libraries(check_fresnel_s ${TR29124_LIBRARIES} quadmath)
run_check(check_fresnel_s)

add_executable(check_gamma_p check/check_gamma_p.cc)
target_include_directories(check_gamma_p PRIVATE .)
target_link_libraries(check_gamma_p ${TR29124_LIBRARIES} quadmath)
run_check(check_gamma_p)

add_executable(check_gamma_q check/check_gamma_q.cc)
target_include_directories(check_gamma_q PRIVATE .)
target_link_libraries(check_gamma_q ${TR29124_LIBRARIES} quadmath)
run_check(check_gamma_q)

add_executable(check_gamma_reciprocal check/check_gamma_reciprocal.cc)
target_include_directories(check_gamma_reciprocal PRIVATE .)
target_link_libraries(check_gamma_reciprocal ${TR29124_LIBRARIES} quadmath)
run_check(check_gamma_reciprocal)

add_executable(check_gegenbauer check/check_gegenbauer.cc)
target_include_directories(check_gegenbauer PRIVATE .)
target_link_libraries(check_gegenbauer ${TR29124_LIBRARIES} quadmath)
run_check(check_gegenbauer)

add_executable(check_hermite check/check_hermite.cc)
target_include_directories(check_hermite PRIVATE .)
target_link_libraries(check_hermite ${TR29124_LIBRARIES} quadmath)
run_check(check_hermite)

add_executable(check_heuman_lambda check/check_heuman_lambda.cc)
target_include_directories(check_heuman_lambda PRIVATE .)
target_link_libraries(check_heuman_lambda ${TR29124_LIBRARIES} quadmath)
run_check(check_heuman_lambda)

add_executable(check_hurwitz_zeta check/check_hurwitz_zeta.cc)
target_include_directories(check_hurwitz_zeta PRIVATE .)
target_link_libraries(check_hurwitz_zeta ${TR29124_LIBRARIES} quadmath)
run_check(check_hurwitz_zeta)

add_executable(check_hyperg check/check_hyperg.cc)
target_include_directories(check_hyperg PRIVATE .)
target_link_libraries(check_hyperg ${TR29124_LIBRARIES} quadmath)
run_check(check_hyperg)

add_executable(check_ibeta check/check_ibeta.cc)
target_include_directories(check_ibeta PRIVATE .)
target_link_libraries(check_ibeta ${TR29124_LIBRARIES} quadmath)
run_check(check_ibeta)

add_executable(check_ibetac check/check_ibetac.cc)
target_include_directories(check_ibetac PRIVATE .)
target_link_libraries(check_ibetac ${TR29124_LIBRARIES} quadmath)
run_check(check_ibetac)

add_executable(check_jacobi check/check_jacobi.cc)
target_include_directories(check_jacobi PRIVATE .)
target_link_libraries(check_jacobi ${TR29124_LIBRARIES} quadmath)

add_executable(check_jacobi_cn check/check_jacobi_cn.cc)
target_include_directories(check_jacobi_cn PRIVATE .)
target_link_libraries(check_jacobi_cn ${TR29124_LIBRARIES} quadmath)
run_check(check_jacobi_cn)

add_executable(check_jacobi_dn check/check_jacobi_dn.cc)
target_include_directories(check_jacobi_dn PRIVATE .)
target_link_libraries(check_jacobi_dn ${TR29124_LIBRARIES} quadmath)
run_check(check_jacobi_dn)

add_executable(check_jacobi_sn check/check_jacobi_sn.cc)
target_include_directories(check_jacobi_sn PRIVATE .)
target_link_libraries(check_jacobi_sn ${TR29124_LIBRARIES} quadmath)
run_check(check_jacobi_sn)

add_executable(check_jacobi_zeta check/check_jacobi_zeta.cc)
target_include_directories(check_jacobi_zeta PRIVATE .)
target_link_libraries(check_jacobi_zeta ${TR29124_LIBRARIES} quadmath)
run_check(check_jacobi_zeta)

add_executable(check_laguerre check/check_laguerre.cc)
target_include_directories(check_laguerre PRIVATE .)
target_link_libraries(check_laguerre ${TR29124_LIBRARIES} quadmath)
run_check(check_laguerre)

add_executable(check_lah check/check_lah.cc)
target_include_directories(check_lah PRIVATE .)
target_link_libraries(check_lah ${TR29124_LIBRARIES} quadmath)
run_check(check_lah)

add_executable(check_lbinomial check/check_lbinomial.cc)
target_include_directories(check_lbinomial PRIVATE .)
target_link_libraries(check_lbinomial ${TR29124_LIBRARIES} quadmath)
run_check(check_lbinomial)

add_executable(check_ldouble_factorial check/check_ldouble_factorial.cc)
target_include_directories(check_ldouble_factorial PRIVATE .)
target_link_libraries(check_ldouble_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_ldouble_factorial)

add_executable(check_legendre check/check_legendre.cc)
target_include_directories(check_legendre PRIVATE .)
target_link_libraries(check_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_legendre)

add_executable(check_legendre_q check/check_legendre_q.cc)
target_include_directories(check_legendre_q PRIVATE .)
target_link_libraries(check_legendre_q ${TR29124_LIBRARIES} quadmath)
run_check(check_legendre_q)

add_executable(check_lfactorial check/check_lfactorial.cc)
target_include_directories(check_lfactorial PRIVATE .)
target_link_libraries(check_lfactorial ${TR29124_LIBRARIES} quadmath)
run_check(check_lfactorial)

add_executable(check_lfalling_factorial check/check_lfalling_factorial.cc)
target_include_directories(check_lfalling_factorial PRIVATE .)
target_link_libraries(check_lfalling_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_lfalling_factorial)

add_executable(check_lgamma check/check_lgamma.cc)
target_include_directories(check_lgamma PRIVATE .)
target_link_libraries(check_lgamma ${TR29124_LIBRARIES} quadmath)
run_check(check_lgamma)

add_executable(check_logistic_p check/check_logistic_p.cc)
target_include_directories(check_logistic_p PRIVATE .)
target_link_libraries(check_logistic_p ${TR29124_LIBRARIES} quadmath)
run_check(check_logistic_p)

add_executable(check_logistic_pdf check/check_logistic_pdf.cc)
target_include_directories(check_logistic_pdf PRIVATE .)
target_link_libraries(check_logistic_pdf ${TR29124_LIBRARIES} quadmath)
run_check(check_logistic_pdf)

add_executable(check_lognormal_p check/check_lognormal_p.cc)
target_include_directories(check_lognormal_p PRIVATE .)
target_link_libraries(check_lognormal_p ${TR29124_LIBRARIES} quadmath)
run_check(check_lognormal_p)

add_executable(check_lognormal_pdf check/check_lognormal_pdf.cc)
target_include_directories(check_lognormal_pdf PRIVATE .)
target_link_libraries(check_lognormal_pdf ${TR29124_LIBRARIES} quadmath)
run_check(check_lognormal_pdf)

add_executable(check_lrising_factorial check/check_lrising_factorial.cc)
target_include_directories(check_lrising_factorial PRIVATE .)
target_link_libraries(check_lrising_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_lrising_factorial)

add_executable(check_normal_p check/check_normal_p.cc)
target_include_directories(check_normal_p PRIVATE .)
target_link_libraries(check_normal_p ${TR29124_LIBRARIES} quadmath)
run_check(check_normal_p)

add_executable(check_normal_pdf check/check_normal_pdf.cc)
target_include_directories(check_normal_pdf PRIVATE .)
target_link_libraries(check_normal_pdf ${TR29124_LIBRARIES} quadmath)
run_check(check_normal_pdf)

add_executable(check_owens_t check/check_owens_t.cc)
target_include_directories(check_owens_t PRIVATE .)
target_link_libraries(check_owens_t ${TR29124_LIBRARIES} quadmath)
run_check(check_owens_t)

add_executable(check_polygamma check/check_polygamma.cc)
target_include_directories(check_polygamma PRIVATE .)
target_link_libraries(check_polygamma ${TR29124_LIBRARIES} quadmath)
run_check(check_polygamma)

add_executable(check_radpoly check/check_radpoly.cc)
target_include_directories(check_radpoly PRIVATE .)
target_link_libraries(check_radpoly ${TR29124_LIBRARIES} quadmath)
run_check(check_radpoly)

add_executable(check_riemann_zeta check/check_riemann_zeta.cc)
target_include_directories(check_riemann_zeta PRIVATE .)
target_link_libraries(check_riemann_zeta ${TR29124_LIBRARIES} quadmath)
run_check(check_riemann_zeta)

add_executable(check_rising_factorial check/check_rising_factorial.cc)
target_include_directories(check_rising_factorial PRIVATE .)
target_link_libraries(check_rising_factorial ${TR29124_LIBRARIES} quadmath)
run_check(check_rising_factorial)

add_executable(check_shi check/check_shi.cc)
target_include_directories(check_shi PRIVATE .)
target_link_libraries(check_shi ${TR29124_LIBRARIES} quadmath)
run_check(check_shi)

add_executable(check_sinc check/check_sinc.cc)
target_include_directories(check_sinc PRIVATE .)
target_link_libraries(check_sinc ${TR29124_LIBRARIES} quadmath)
run_check(check_sinc)

add_executable(check_sinc_pi check/check_sinc_pi.cc)
target_include_directories(check_sinc_pi PRIVATE .)
target_link_libraries(check_sinc_pi ${TR29124_LIBRARIES} quadmath)
run_check(check_sinc_pi)

add_executable(check_sinhint check/check_sinhint.cc)
target_include_directories(check_sinhint PRIVATE .)
target_link_libraries(check_sinhint ${TR29124_LIBRARIES} quadmath)
run_check(check_sinhint)

add_executable(check_sinint check/check_sinint.cc)
target_include_directories(check_sinint PRIVATE .)
target_link_libraries(check_sinint ${TR29124_LIBRARIES} quadmath)
run_check(check_sinint)

add_executable(check_sin_pi check/check_sin_pi.cc)
target_include_directories(check_sin_pi PRIVATE .)
target_link_libraries(check_sin_pi ${TR29124_LIBRARIES} quadmath)
run_check(check_sin_pi)

add_executable(check_sph_bessel check/check_sph_bessel.cc)
target_include_directories(check_sph_bessel PRIVATE .)
target_link_libraries(check_sph_bessel ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_bessel)

add_executable(check_sph_bessel_i check/check_sph_bessel_i.cc)
target_include_directories(check_sph_bessel_i PRIVATE .)
target_link_libraries(check_sph_bessel_i ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_bessel_i)

add_executable(check_sph_bessel_k check/check_sph_bessel_k.cc)
target_include_directories(check_sph_bessel_k PRIVATE .)
target_link_libraries(check_sph_bessel_k ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_bessel_k)

add_executable(check_sph_hankel_1 check/check_sph_hankel_1.cc)
target_include_directories(check_sph_hankel_1 PRIVATE .)
target_link_libraries(check_sph_hankel_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_hankel_1)

add_executable(check_sph_hankel_2 check/check_sph_hankel_2.cc)
target_include_directories(check_sph_hankel_2 PRIVATE .)
target_link_libraries(check_sph_hankel_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_hankel_2)

add_executable(check_sph_harmonic check/check_sph_harmonic.cc)
target_include_directories(check_sph_harmonic PRIVATE .)
target_link_libraries(check_sph_harmonic ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_harmonic)

add_executable(check_sph_legendre check/check_sph_legendre.cc)
target_include_directories(check_sph_legendre PRIVATE .)
target_link_libraries(check_sph_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_legendre)

add_executable(check_sph_neumann check/check_sph_neumann.cc)
target_include_directories(check_sph_neumann PRIVATE .)
target_link_libraries(check_sph_neumann ${TR29124_LIBRARIES} quadmath)
run_check(check_sph_neumann)

add_executable(check_stirling_1 check/check_stirling_1.cc)
target_include_directories(check_stirling_1 PRIVATE .)
target_link_libraries(check_stirling_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_stirling_1)

add_executable(check_stirling_2 check/check_stirling_2.cc)
target_include_directories(check_stirling_2 PRIVATE .)
target_link_libraries(check_stirling_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_stirling_2)

add_executable(check_tgamma check/check_tgamma.cc)
target_include_directories(check_tgamma PRIVATE .)
target_link_libraries(check_tgamma ${TR29124_LIBRARIES} quadmath)
run_check(check_tgamma)

add_executable(check_tgamma_lower check/check_tgamma_lower.cc)
target_include_directories(check_tgamma_lower PRIVATE .)
target_link_libraries(check_tgamma_lower ${TR29124_LIBRARIES} quadmath)
run_check(check_tgamma_lower)

add_executable(check_theta_1 check/check_theta_1.cc)
target_include_directories(check_theta_1 PRIVATE .)
target_link_libraries(check_theta_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_theta_1)

add_executable(check_theta_2 check/check_theta_2.cc)
target_include_directories(check_theta_2 PRIVATE .)
target_link_libraries(check_theta_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_theta_2)

add_executable(check_theta_3 check/check_theta_3.cc)
target_include_directories(check_theta_3 PRIVATE .)
target_link_libraries(check_theta_3 ${TR29124_LIBRARIES} quadmath)
run_check(check_theta_3)

add_executable(check_theta_4 check/check_theta_4.cc)
target_include_directories(check_theta_4 PRIVATE .)
target_link_libraries(check_theta_4 ${TR29124_LIBRARIES} quadmath)
run_check(check_theta_4)

add_executable(check_tr1_assoc_laguerre check/check_tr1_assoc_laguerre.cc)
target_include_directories(check_tr1_assoc_laguerre PRIVATE .)
target_link_libraries(check_tr1_assoc_laguerre ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_assoc_laguerre)

add_executable(check_tr1_assoc_legendre check/check_tr1_assoc_legendre.cc)
target_include_directories(check_tr1_assoc_legendre PRIVATE .)
target_link_libraries(check_tr1_assoc_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_assoc_legendre)

add_executable(check_tr1_beta check/check_tr1_beta.cc)
target_include_directories(check_tr1_beta PRIVATE .)
target_link_libraries(check_tr1_beta ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_beta)

add_executable(check_tr1_comp_ellint_1 check/check_tr1_comp_ellint_1.cc)
target_include_directories(check_tr1_comp_ellint_1 PRIVATE .)
target_link_libraries(check_tr1_comp_ellint_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_comp_ellint_1)

add_executable(check_tr1_comp_ellint_2 check/check_tr1_comp_ellint_2.cc)
target_include_directories(check_tr1_comp_ellint_2 PRIVATE .)
target_link_libraries(check_tr1_comp_ellint_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_comp_ellint_2)

add_executable(check_tr1_comp_ellint_3 check/check_tr1_comp_ellint_3.cc)
target_include_directories(check_tr1_comp_ellint_3 PRIVATE .)
target_link_libraries(check_tr1_comp_ellint_3 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_comp_ellint_3)

add_executable(check_tr1_conf_hyperg check/check_tr1_conf_hyperg.cc)
target_include_directories(check_tr1_conf_hyperg PRIVATE .)
target_link_libraries(check_tr1_conf_hyperg ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_conf_hyperg)

add_executable(check_tr1_cyl_bessel_i check/check_tr1_cyl_bessel_i.cc)
target_include_directories(check_tr1_cyl_bessel_i PRIVATE .)
target_link_libraries(check_tr1_cyl_bessel_i ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_cyl_bessel_i)

add_executable(check_tr1_cyl_bessel_j check/check_tr1_cyl_bessel_j.cc)
target_include_directories(check_tr1_cyl_bessel_j PRIVATE .)
target_link_libraries(check_tr1_cyl_bessel_j ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_cyl_bessel_j)

add_executable(check_tr1_cyl_bessel_k check/check_tr1_cyl_bessel_k.cc)
target_include_directories(check_tr1_cyl_bessel_k PRIVATE .)
target_link_libraries(check_tr1_cyl_bessel_k ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_cyl_bessel_k)

add_executable(check_tr1_cyl_neumann check/check_tr1_cyl_neumann.cc)
target_include_directories(check_tr1_cyl_neumann PRIVATE .)
target_link_libraries(check_tr1_cyl_neumann ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_cyl_neumann)

add_executable(check_tr1_ellint_1 check/check_tr1_ellint_1.cc)
target_include_directories(check_tr1_ellint_1 PRIVATE .)
target_link_libraries(check_tr1_ellint_1 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_ellint_1)

add_executable(check_tr1_ellint_2 check/check_tr1_ellint_2.cc)
target_include_directories(check_tr1_ellint_2 PRIVATE .)
target_link_libraries(check_tr1_ellint_2 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_ellint_2)

add_executable(check_tr1_ellint_3 check/check_tr1_ellint_3.cc)
target_include_directories(check_tr1_ellint_3 PRIVATE .)
target_link_libraries(check_tr1_ellint_3 ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_ellint_3)

add_executable(check_tr1_expint check/check_tr1_expint.cc)
target_include_directories(check_tr1_expint PRIVATE .)
target_link_libraries(check_tr1_expint ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_expint)

add_executable(check_tr1_hermite check/check_tr1_hermite.cc)
target_include_directories(check_tr1_hermite PRIVATE .)
target_link_libraries(check_tr1_hermite ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_hermite)

add_executable(check_tr1_hyperg check/check_tr1_hyperg.cc)
target_include_directories(check_tr1_hyperg PRIVATE .)
target_link_libraries(check_tr1_hyperg ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_hyperg)

add_executable(check_tr1_laguerre check/check_tr1_laguerre.cc)
target_include_directories(check_tr1_laguerre PRIVATE .)
target_link_libraries(check_tr1_laguerre ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_laguerre)

add_executable(check_tr1_legendre check/check_tr1_legendre.cc)
target_include_directories(check_tr1_legendre PRIVATE .)
target_link_libraries(check_tr1_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_legendre)

add_executable(check_tr1_riemann_zeta check/check_tr1_riemann_zeta.cc)
target_include_directories(check_tr1_riemann_zeta PRIVATE .)
target_link_libraries(check_tr1_riemann_zeta ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_riemann_zeta)

add_executable(check_tr1_sph_bessel check/check_tr1_sph_bessel.cc)
target_include_directories(check_tr1_sph_bessel PRIVATE .)
target_link_libraries(check_tr1_sph_bessel ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_sph_bessel)

add_executable(check_tr1_sph_legendre check/check_tr1_sph_legendre.cc)
target_include_directories(check_tr1_sph_legendre PRIVATE .)
target_link_libraries(check_tr1_sph_legendre ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_sph_legendre)

add_executable(check_tr1_sph_neumann check/check_tr1_sph_neumann.cc)
target_include_directories(check_tr1_sph_neumann PRIVATE .)
target_link_libraries(check_tr1_sph_neumann ${TR29124_LIBRARIES} quadmath)
run_check(check_tr1_sph_neumann)

add_executable(check_zernike check/check_zernike.cc)
target_include_directories(check_zernike PRIVATE .)
target_link_libraries(check_zernike ${TR29124_LIBRARIES} quadmath)
run_check(check_zernike)

add_executable(complex_airy_ai check/complex_airy_ai.cc)
target_include_directories(complex_airy_ai PRIVATE .)
target_link_libraries(complex_airy_ai ${TR29124_LIBRARIES} quadmath)
run_check(complex_airy_ai)

add_executable(complex_airy_bi check/complex_airy_bi.cc)
target_include_directories(complex_airy_bi PRIVATE .)
target_link_libraries(complex_airy_bi ${TR29124_LIBRARIES} quadmath)
run_check(complex_airy_bi)

add_executable(complex_ellint_rc check/complex_ellint_rc.cc)
target_include_directories(complex_ellint_rc PRIVATE .)
target_link_libraries(complex_ellint_rc ${TR29124_LIBRARIES} quadmath)
run_check(complex_ellint_rc)

add_executable(complex_ellint_rd check/complex_ellint_rd.cc)
target_include_directories(complex_ellint_rd PRIVATE .)
target_link_libraries(complex_ellint_rd ${TR29124_LIBRARIES} quadmath)
run_check(complex_ellint_rd)

add_executable(complex_ellint_rf check/complex_ellint_rf.cc)
target_include_directories(complex_ellint_rf PRIVATE .)
target_link_libraries(complex_ellint_rf ${TR29124_LIBRARIES} quadmath)
run_check(complex_ellint_rf)

add_executable(complex_ellint_rg check/complex_ellint_rg.cc)
target_include_directories(complex_ellint_rg PRIVATE .)
target_link_libraries(complex_ellint_rg ${TR29124_LIBRARIES} quadmath)
run_check(complex_ellint_rg)

add_executable(complex_ellint_rj check/complex_ellint_rj.cc)
target_include_directories(complex_ellint_rj PRIVATE .)
target_link_libraries(complex_ellint_rj ${TR29124_LIBRARIES} quadmath)
run_check(complex_ellint_rj)

add_executable(deathmatch_comp_ellint check/deathmatch_comp_ellint.cc)
target_include_directories(deathmatch_comp_ellint PRIVATE .)
target_link_libraries(deathmatch_comp_ellint ${TR29124_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg check/deathmatch_conf_hyperg.cc)
target_include_directories(deathmatch_conf_hyperg PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg ${TR29124_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg_lim check/deathmatch_conf_hyperg_lim.cc)
target_include_directories(deathmatch_conf_hyperg_lim PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg_lim ${TR29124_LIBRARIES} quadmath)
run_check(deathmatch_conf_hyperg_lim)

add_executable(deathmatch_hyperg check/deathmatch_hyperg.cc)
target_include_directories(deathmatch_hyperg PRIVATE .)
target_link_libraries(deathmatch_hyperg ${TR29124_LIBRARIES} quadmath)
run_check(deathmatch_hyperg)

add_executable(origin_cyl_bessel_j check/origin_cyl_bessel_j.cc)
target_include_directories(origin_cyl_bessel_j PRIVATE .)
target_link_libraries(origin_cyl_bessel_j ${TR29124_LIBRARIES} quadmath)
run_check(origin_cyl_bessel_j)

add_executable(origin_cyl_neumann check/origin_cyl_neumann.cc)
target_include_directories(origin_cyl_neumann PRIVATE .)
target_link_libraries(origin_cyl_neumann ${TR29124_LIBRARIES} quadmath)
run_check(origin_cyl_neumann)

add_executable(pr56216_cyl_bessel_i check/pr56216_cyl_bessel_i.cc)
target_include_directories(pr56216_cyl_bessel_i PRIVATE .)
target_link_libraries(pr56216_cyl_bessel_i ${TR29124_LIBRARIES} quadmath)
run_check(pr56216_cyl_bessel_i)

add_executable(pr56216_cyl_hankel_1 check/pr56216_cyl_hankel_1.cc)
target_include_directories(pr56216_cyl_hankel_1 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_1 ${TR29124_LIBRARIES} quadmath)
run_check(pr56216_cyl_hankel_1)

add_executable(pr56216_cyl_hankel_2 check/pr56216_cyl_hankel_2.cc)
target_include_directories(pr56216_cyl_hankel_2 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_2 ${TR29124_LIBRARIES} quadmath)
run_check(pr56216_cyl_hankel_2)

add_executable(pr68397 check/pr68397.cc)
target_include_directories(pr68397 PRIVATE .)
target_link_libraries(pr68397 ${TR29124_LIBRARIES} quadmath)
run_check(pr68397)

add_executable(pr83140 check/pr83140.cc)
target_include_directories(pr83140 PRIVATE .)
target_link_libraries(pr83140 ${TR29124_LIBRARIES} quadmath)
run_check(pr83140)

add_executable(pr86655_assoc_legendre check/pr86655_assoc_legendre.cc)
target_include_directories(pr86655_assoc_legendre PRIVATE .)
target_link_libraries(pr86655_assoc_legendre ${TR29124_LIBRARIES} quadmath)
run_check(pr86655_assoc_legendre)

add_executable(pr86655_sph_legendre check/pr86655_sph_legendre.cc)
target_include_directories(pr86655_sph_legendre PRIVATE .)
target_link_libraries(pr86655_sph_legendre ${TR29124_LIBRARIES} quadmath)
run_check(pr86655_sph_legendre)

add_executable(pr86655_tr1_assoc_legendre check/pr86655_tr1_assoc_legendre.cc)
target_include_directories(pr86655_tr1_assoc_legendre PRIVATE .)
target_link_libraries(pr86655_tr1_assoc_legendre ${TR29124_LIBRARIES} quadmath)
run_check(pr86655_tr1_assoc_legendre)

add_executable(pr86655_tr1_sph_legendre check/pr86655_tr1_sph_legendre.cc)
target_include_directories(pr86655_tr1_sph_legendre PRIVATE .)
target_link_libraries(pr86655_tr1_sph_legendre ${TR29124_LIBRARIES} quadmath)
run_check(pr86655_tr1_sph_legendre)

# This needs C++14!
#add_executable(tr1_airy check/tr1_airy.cc)
#target_include_directories(tr1_airy PRIVATE .)
#target_link_libraries(tr1_airy ${TR29124_LIBRARIES} quadmath)
#run_check(tr1_airy)

add_executable(tr1_pr68397 check/tr1_pr68397.cc)
target_include_directories(tr1_pr68397 PRIVATE .)
target_link_libraries(tr1_pr68397 ${TR29124_LIBRARIES} quadmath)
run_check(tr1_pr68397)

add_executable(tr1_pr83140 check/tr1_pr83140.cc)
target_include_directories(tr1_pr83140 PRIVATE .)
target_link_libraries(tr1_pr83140 ${TR29124_LIBRARIES} quadmath)
run_check(tr1_pr83140)

