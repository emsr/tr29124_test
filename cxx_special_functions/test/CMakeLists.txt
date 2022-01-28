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

add_custom_target(make_libstdcpp_support_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_definitions(-D__TEST_DEBUG)

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
  cxx_integration)

add_executable(testcase testcase.cpp)
target_link_libraries(testcase wrap_gsl wrap_boost wrap_burkhardt)
target_include_directories(testcase PRIVATE .)

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

add_executable(check_airy_ai src/check_airy_ai.cc)
target_include_directories(check_airy_ai PRIVATE .)
target_link_libraries(check_airy_ai ${EMSR_LIBRARIES} quadmath)
run_check(check_airy_ai)

add_executable(check_airy_bi src/check_airy_bi.cc)
target_include_directories(check_airy_bi PRIVATE .)
target_link_libraries(check_airy_bi ${EMSR_LIBRARIES} quadmath)
run_check(check_airy_bi)

add_executable(check_assoc_laguerre src/check_assoc_laguerre.cc)
target_include_directories(check_assoc_laguerre PRIVATE .)
target_link_libraries(check_assoc_laguerre ${EMSR_LIBRARIES} quadmath)
run_check(check_assoc_laguerre)

add_executable(check_assoc_legendre src/check_assoc_legendre.cc)
target_include_directories(check_assoc_legendre PRIVATE .)
target_link_libraries(check_assoc_legendre ${EMSR_LIBRARIES} quadmath)
run_check(check_assoc_legendre)

add_executable(check_bell src/check_bell.cc)
target_include_directories(check_bell PRIVATE .)
target_link_libraries(check_bell ${EMSR_LIBRARIES} quadmath)
run_check(check_bell)

add_executable(check_bernoulli src/check_bernoulli.cc)
target_include_directories(check_bernoulli PRIVATE .)
target_link_libraries(check_bernoulli ${EMSR_LIBRARIES} quadmath)
run_check(check_bernoulli)

add_executable(check_beta src/check_beta.cc)
target_include_directories(check_beta PRIVATE .)
target_link_libraries(check_beta ${EMSR_LIBRARIES} quadmath)
run_check(check_beta)

add_executable(check_binomial src/check_binomial.cc)
target_include_directories(check_binomial PRIVATE .)
target_link_libraries(check_binomial ${EMSR_LIBRARIES} quadmath)
run_check(check_binomial)

add_executable(check_chebyshev_t src/check_chebyshev_t.cc)
target_include_directories(check_chebyshev_t PRIVATE .)
target_link_libraries(check_chebyshev_t ${EMSR_LIBRARIES} quadmath)
run_check(check_chebyshev_t)

add_executable(check_chebyshev_u src/check_chebyshev_u.cc)
target_include_directories(check_chebyshev_u PRIVATE .)
target_link_libraries(check_chebyshev_u ${EMSR_LIBRARIES} quadmath)
run_check(check_chebyshev_u)

add_executable(check_chebyshev_v src/check_chebyshev_v.cc)
target_include_directories(check_chebyshev_v PRIVATE .)
target_link_libraries(check_chebyshev_v ${EMSR_LIBRARIES} quadmath)
run_check(check_chebyshev_v)

add_executable(check_chebyshev_w src/check_chebyshev_w.cc)
target_include_directories(check_chebyshev_w PRIVATE .)
target_link_libraries(check_chebyshev_w ${EMSR_LIBRARIES} quadmath)
run_check(check_chebyshev_w)

add_executable(check_chi src/check_chi.cc)
target_include_directories(check_chi PRIVATE .)
target_link_libraries(check_chi ${EMSR_LIBRARIES} quadmath)
run_check(check_chi)

add_executable(check_clausen_cl src/check_clausen_cl.cc)
target_include_directories(check_clausen_cl PRIVATE .)
target_link_libraries(check_clausen_cl ${EMSR_LIBRARIES} quadmath)
run_check(check_clausen_cl)

add_executable(check_comp_ellint_1 src/check_comp_ellint_1.cc)
target_include_directories(check_comp_ellint_1 PRIVATE .)
target_link_libraries(check_comp_ellint_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_comp_ellint_1)

add_executable(check_comp_ellint_2 src/check_comp_ellint_2.cc)
target_include_directories(check_comp_ellint_2 PRIVATE .)
target_link_libraries(check_comp_ellint_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_comp_ellint_2)

add_executable(check_comp_ellint_3 src/check_comp_ellint_3.cc)
target_include_directories(check_comp_ellint_3 PRIVATE .)
target_link_libraries(check_comp_ellint_3 ${EMSR_LIBRARIES} quadmath)
run_check(check_comp_ellint_3)

add_executable(check_comp_ellint_d src/check_comp_ellint_d.cc)
target_include_directories(check_comp_ellint_d PRIVATE .)
target_link_libraries(check_comp_ellint_d ${EMSR_LIBRARIES} quadmath)
run_check(check_comp_ellint_d)

add_executable(check_conf_hyperg src/check_conf_hyperg.cc)
target_include_directories(check_conf_hyperg PRIVATE .)
target_link_libraries(check_conf_hyperg ${EMSR_LIBRARIES} quadmath)
run_check(check_conf_hyperg)

add_executable(check_conf_hyperg_lim src/check_conf_hyperg_lim.cc)
target_include_directories(check_conf_hyperg_lim PRIVATE .)
target_link_libraries(check_conf_hyperg_lim ${EMSR_LIBRARIES} quadmath)
run_check(check_conf_hyperg_lim)

add_executable(check_coshint src/check_coshint.cc)
target_include_directories(check_coshint PRIVATE .)
target_link_libraries(check_coshint ${EMSR_LIBRARIES} quadmath)
run_check(check_coshint)

add_executable(check_cosint src/check_cosint.cc)
target_include_directories(check_cosint PRIVATE .)
target_link_libraries(check_cosint ${EMSR_LIBRARIES} quadmath)
run_check(check_cosint)

add_executable(check_cos_pi src/check_cos_pi.cc)
target_include_directories(check_cos_pi PRIVATE .)
target_link_libraries(check_cos_pi ${EMSR_LIBRARIES} quadmath)
run_check(check_cos_pi)

add_executable(check_cyl_bessel_i src/check_cyl_bessel_i.cc)
target_include_directories(check_cyl_bessel_i PRIVATE .)
target_link_libraries(check_cyl_bessel_i ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_bessel_i)

add_executable(check_cyl_bessel_i_scaled src/check_cyl_bessel_i_scaled.cc)
target_include_directories(check_cyl_bessel_i_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_i_scaled ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_bessel_i_scaled)

add_executable(check_cyl_bessel_j src/check_cyl_bessel_j.cc)
target_include_directories(check_cyl_bessel_j PRIVATE .)
target_link_libraries(check_cyl_bessel_j ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_bessel_j)

add_executable(check_cyl_bessel_k src/check_cyl_bessel_k.cc)
target_include_directories(check_cyl_bessel_k PRIVATE .)
target_link_libraries(check_cyl_bessel_k ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_bessel_k)

add_executable(check_cyl_bessel_k_scaled src/check_cyl_bessel_k_scaled.cc)
target_include_directories(check_cyl_bessel_k_scaled PRIVATE .)
target_link_libraries(check_cyl_bessel_k_scaled ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_bessel_k_scaled)

add_executable(check_cyl_hankel_1 src/check_cyl_hankel_1.cc)
target_include_directories(check_cyl_hankel_1 PRIVATE .)
target_link_libraries(check_cyl_hankel_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_hankel_1)

add_executable(check_cyl_hankel_2 src/check_cyl_hankel_2.cc)
target_include_directories(check_cyl_hankel_2 PRIVATE .)
target_link_libraries(check_cyl_hankel_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_hankel_2)

add_executable(check_cyl_neumann src/check_cyl_neumann.cc)
target_include_directories(check_cyl_neumann PRIVATE .)
target_link_libraries(check_cyl_neumann ${EMSR_LIBRARIES} quadmath)
run_check(check_cyl_neumann)

add_executable(check_dawson src/check_dawson.cc)
target_include_directories(check_dawson PRIVATE .)
target_link_libraries(check_dawson ${EMSR_LIBRARIES} quadmath)
run_check(check_dawson)

add_executable(check_debye src/check_debye.cc)
target_include_directories(check_debye PRIVATE .)
target_link_libraries(check_debye ${EMSR_LIBRARIES} quadmath)
run_check(check_debye)

add_executable(check_digamma src/check_digamma.cc)
target_include_directories(check_digamma PRIVATE .)
target_link_libraries(check_digamma ${EMSR_LIBRARIES} quadmath)
run_check(check_digamma)

add_executable(check_dilog src/check_dilog.cc)
target_include_directories(check_dilog PRIVATE .)
target_link_libraries(check_dilog ${EMSR_LIBRARIES} quadmath)
run_check(check_dilog)

add_executable(check_dirichlet_beta src/check_dirichlet_beta.cc)
target_include_directories(check_dirichlet_beta PRIVATE .)
target_link_libraries(check_dirichlet_beta ${EMSR_LIBRARIES} quadmath)
run_check(check_dirichlet_beta)

add_executable(check_dirichlet_eta src/check_dirichlet_eta.cc)
target_include_directories(check_dirichlet_eta PRIVATE .)
target_link_libraries(check_dirichlet_eta ${EMSR_LIBRARIES} quadmath)
run_check(check_dirichlet_eta)

add_executable(check_dirichlet_lambda src/check_dirichlet_lambda.cc)
target_include_directories(check_dirichlet_lambda PRIVATE .)
target_link_libraries(check_dirichlet_lambda ${EMSR_LIBRARIES} quadmath)
run_check(check_dirichlet_lambda)

add_executable(check_double_factorial src/check_double_factorial.cc)
target_include_directories(check_double_factorial PRIVATE .)
target_link_libraries(check_double_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_double_factorial)

add_executable(check_ellint_1 src/check_ellint_1.cc)
target_include_directories(check_ellint_1 PRIVATE .)
target_link_libraries(check_ellint_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_1)

add_executable(check_ellint_2 src/check_ellint_2.cc)
target_include_directories(check_ellint_2 PRIVATE .)
target_link_libraries(check_ellint_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_2)

add_executable(check_ellint_3 src/check_ellint_3.cc)
target_include_directories(check_ellint_3 PRIVATE .)
target_link_libraries(check_ellint_3 ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_3)

add_executable(check_ellint_d src/check_ellint_d.cc)
target_include_directories(check_ellint_d PRIVATE .)
target_link_libraries(check_ellint_d ${EMSR_LIBRARIES} quadmath)

add_executable(check_ellint_rc src/check_ellint_rc.cc)
target_include_directories(check_ellint_rc PRIVATE .)
target_link_libraries(check_ellint_rc ${EMSR_LIBRARIES} quadmath)

add_executable(check_ellint_rd src/check_ellint_rd.cc)
target_include_directories(check_ellint_rd PRIVATE .)
target_link_libraries(check_ellint_rd ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_rd)

add_executable(check_ellint_rf src/check_ellint_rf.cc)
target_include_directories(check_ellint_rf PRIVATE .)
target_link_libraries(check_ellint_rf ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_rf)

add_executable(check_ellint_rg src/check_ellint_rg.cc)
target_include_directories(check_ellint_rg PRIVATE .)
target_link_libraries(check_ellint_rg ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_rg)

add_executable(check_ellint_rj src/check_ellint_rj.cc)
target_include_directories(check_ellint_rj PRIVATE .)
target_link_libraries(check_ellint_rj ${EMSR_LIBRARIES} quadmath)
run_check(check_ellint_rj)

add_executable(check_ellnome src/check_ellnome.cc)
target_include_directories(check_ellnome PRIVATE .)
target_link_libraries(check_ellnome ${EMSR_LIBRARIES} quadmath)
run_check(check_ellnome)

add_executable(check_euler src/check_euler.cc)
target_include_directories(check_euler PRIVATE .)
target_link_libraries(check_euler ${EMSR_LIBRARIES} quadmath)
run_check(check_euler)

add_executable(check_eulerian_1 src/check_eulerian_1.cc)
target_include_directories(check_eulerian_1 PRIVATE .)
target_link_libraries(check_eulerian_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_eulerian_1)

add_executable(check_eulerian_2 src/check_eulerian_2.cc)
target_include_directories(check_eulerian_2 PRIVATE .)
target_link_libraries(check_eulerian_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_eulerian_2)

add_executable(check_expint src/check_expint.cc)
target_include_directories(check_expint PRIVATE .)
target_link_libraries(check_expint ${EMSR_LIBRARIES} quadmath)
run_check(check_expint)

add_executable(check_expint_en src/check_expint_en.cc)
target_include_directories(check_expint_en PRIVATE .)
target_link_libraries(check_expint_en ${EMSR_LIBRARIES} quadmath)
run_check(check_expint_en)

add_executable(check_factorial src/check_factorial.cc)
target_include_directories(check_factorial PRIVATE .)
target_link_libraries(check_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_factorial)

add_executable(check_falling_factorial src/check_falling_factorial.cc)
target_include_directories(check_falling_factorial PRIVATE .)
target_link_libraries(check_falling_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_falling_factorial)

add_executable(check_fresnel_c src/check_fresnel_c.cc)
target_include_directories(check_fresnel_c PRIVATE .)
target_link_libraries(check_fresnel_c ${EMSR_LIBRARIES} quadmath)
run_check(check_fresnel_c)

add_executable(check_fresnel_s src/check_fresnel_s.cc)
target_include_directories(check_fresnel_s PRIVATE .)
target_link_libraries(check_fresnel_s ${EMSR_LIBRARIES} quadmath)
run_check(check_fresnel_s)

add_executable(check_gamma_p src/check_gamma_p.cc)
target_include_directories(check_gamma_p PRIVATE .)
target_link_libraries(check_gamma_p ${EMSR_LIBRARIES} quadmath)
run_check(check_gamma_p)

add_executable(check_gamma_q src/check_gamma_q.cc)
target_include_directories(check_gamma_q PRIVATE .)
target_link_libraries(check_gamma_q ${EMSR_LIBRARIES} quadmath)
run_check(check_gamma_q)

add_executable(check_gamma_reciprocal src/check_gamma_reciprocal.cc)
target_include_directories(check_gamma_reciprocal PRIVATE .)
target_link_libraries(check_gamma_reciprocal ${EMSR_LIBRARIES} quadmath)
run_check(check_gamma_reciprocal)

add_executable(check_gegenbauer src/check_gegenbauer.cc)
target_include_directories(check_gegenbauer PRIVATE .)
target_link_libraries(check_gegenbauer ${EMSR_LIBRARIES} quadmath)
run_check(check_gegenbauer)

add_executable(check_hermite src/check_hermite.cc)
target_include_directories(check_hermite PRIVATE .)
target_link_libraries(check_hermite ${EMSR_LIBRARIES} quadmath)
run_check(check_hermite)

add_executable(check_heuman_lambda src/check_heuman_lambda.cc)
target_include_directories(check_heuman_lambda PRIVATE .)
target_link_libraries(check_heuman_lambda ${EMSR_LIBRARIES} quadmath)
run_check(check_heuman_lambda)

add_executable(check_hurwitz_zeta src/check_hurwitz_zeta.cc)
target_include_directories(check_hurwitz_zeta PRIVATE .)
target_link_libraries(check_hurwitz_zeta ${EMSR_LIBRARIES} quadmath)
run_check(check_hurwitz_zeta)

add_executable(check_hyperg src/check_hyperg.cc)
target_include_directories(check_hyperg PRIVATE .)
target_link_libraries(check_hyperg ${EMSR_LIBRARIES} quadmath)
run_check(check_hyperg)

add_executable(check_ibeta src/check_ibeta.cc)
target_include_directories(check_ibeta PRIVATE .)
target_link_libraries(check_ibeta ${EMSR_LIBRARIES} quadmath)
run_check(check_ibeta)

add_executable(check_ibetac src/check_ibetac.cc)
target_include_directories(check_ibetac PRIVATE .)
target_link_libraries(check_ibetac ${EMSR_LIBRARIES} quadmath)
run_check(check_ibetac)

add_executable(check_jacobi src/check_jacobi.cc)
target_include_directories(check_jacobi PRIVATE .)
target_link_libraries(check_jacobi ${EMSR_LIBRARIES} quadmath)

add_executable(check_jacobi_cn src/check_jacobi_cn.cc)
target_include_directories(check_jacobi_cn PRIVATE .)
target_link_libraries(check_jacobi_cn ${EMSR_LIBRARIES} quadmath)
run_check(check_jacobi_cn)

add_executable(check_jacobi_dn src/check_jacobi_dn.cc)
target_include_directories(check_jacobi_dn PRIVATE .)
target_link_libraries(check_jacobi_dn ${EMSR_LIBRARIES} quadmath)
run_check(check_jacobi_dn)

add_executable(check_jacobi_sn src/check_jacobi_sn.cc)
target_include_directories(check_jacobi_sn PRIVATE .)
target_link_libraries(check_jacobi_sn ${EMSR_LIBRARIES} quadmath)
run_check(check_jacobi_sn)

add_executable(check_jacobi_zeta src/check_jacobi_zeta.cc)
target_include_directories(check_jacobi_zeta PRIVATE .)
target_link_libraries(check_jacobi_zeta ${EMSR_LIBRARIES} quadmath)
run_check(check_jacobi_zeta)

add_executable(check_laguerre src/check_laguerre.cc)
target_include_directories(check_laguerre PRIVATE .)
target_link_libraries(check_laguerre ${EMSR_LIBRARIES} quadmath)
run_check(check_laguerre)

add_executable(check_lah src/check_lah.cc)
target_include_directories(check_lah PRIVATE .)
target_link_libraries(check_lah ${EMSR_LIBRARIES} quadmath)
run_check(check_lah)

add_executable(check_lbinomial src/check_lbinomial.cc)
target_include_directories(check_lbinomial PRIVATE .)
target_link_libraries(check_lbinomial ${EMSR_LIBRARIES} quadmath)
run_check(check_lbinomial)

add_executable(check_ldouble_factorial src/check_ldouble_factorial.cc)
target_include_directories(check_ldouble_factorial PRIVATE .)
target_link_libraries(check_ldouble_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_ldouble_factorial)

add_executable(check_legendre src/check_legendre.cc)
target_include_directories(check_legendre PRIVATE .)
target_link_libraries(check_legendre ${EMSR_LIBRARIES} quadmath)
run_check(check_legendre)

add_executable(check_legendre_q src/check_legendre_q.cc)
target_include_directories(check_legendre_q PRIVATE .)
target_link_libraries(check_legendre_q ${EMSR_LIBRARIES} quadmath)
run_check(check_legendre_q)

add_executable(check_lfactorial src/check_lfactorial.cc)
target_include_directories(check_lfactorial PRIVATE .)
target_link_libraries(check_lfactorial ${EMSR_LIBRARIES} quadmath)
run_check(check_lfactorial)

add_executable(check_lfalling_factorial src/check_lfalling_factorial.cc)
target_include_directories(check_lfalling_factorial PRIVATE .)
target_link_libraries(check_lfalling_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_lfalling_factorial)

add_executable(check_lgamma src/check_lgamma.cc)
target_include_directories(check_lgamma PRIVATE .)
target_link_libraries(check_lgamma ${EMSR_LIBRARIES} quadmath)
run_check(check_lgamma)

add_executable(check_logistic_p src/check_logistic_p.cc)
target_include_directories(check_logistic_p PRIVATE .)
target_link_libraries(check_logistic_p ${EMSR_LIBRARIES} quadmath)
run_check(check_logistic_p)

add_executable(check_logistic_pdf src/check_logistic_pdf.cc)
target_include_directories(check_logistic_pdf PRIVATE .)
target_link_libraries(check_logistic_pdf ${EMSR_LIBRARIES} quadmath)
run_check(check_logistic_pdf)

add_executable(check_lognormal_p src/check_lognormal_p.cc)
target_include_directories(check_lognormal_p PRIVATE .)
target_link_libraries(check_lognormal_p ${EMSR_LIBRARIES} quadmath)
run_check(check_lognormal_p)

add_executable(check_lognormal_pdf src/check_lognormal_pdf.cc)
target_include_directories(check_lognormal_pdf PRIVATE .)
target_link_libraries(check_lognormal_pdf ${EMSR_LIBRARIES} quadmath)
run_check(check_lognormal_pdf)

add_executable(check_lrising_factorial src/check_lrising_factorial.cc)
target_include_directories(check_lrising_factorial PRIVATE .)
target_link_libraries(check_lrising_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_lrising_factorial)

add_executable(check_normal_p src/check_normal_p.cc)
target_include_directories(check_normal_p PRIVATE .)
target_link_libraries(check_normal_p ${EMSR_LIBRARIES} quadmath)
run_check(check_normal_p)

add_executable(check_normal_pdf src/check_normal_pdf.cc)
target_include_directories(check_normal_pdf PRIVATE .)
target_link_libraries(check_normal_pdf ${EMSR_LIBRARIES} quadmath)
run_check(check_normal_pdf)

add_executable(check_owens_t src/check_owens_t.cc)
target_include_directories(check_owens_t PRIVATE .)
target_link_libraries(check_owens_t ${EMSR_LIBRARIES} quadmath)
run_check(check_owens_t)

add_executable(check_polygamma src/check_polygamma.cc)
target_include_directories(check_polygamma PRIVATE .)
target_link_libraries(check_polygamma ${EMSR_LIBRARIES} quadmath)
run_check(check_polygamma)

add_executable(check_radpoly src/check_radpoly.cc)
target_include_directories(check_radpoly PRIVATE .)
target_link_libraries(check_radpoly ${EMSR_LIBRARIES} quadmath)
run_check(check_radpoly)

add_executable(check_riemann_zeta src/check_riemann_zeta.cc)
target_include_directories(check_riemann_zeta PRIVATE .)
target_link_libraries(check_riemann_zeta ${EMSR_LIBRARIES} quadmath)
run_check(check_riemann_zeta)

add_executable(check_rising_factorial src/check_rising_factorial.cc)
target_include_directories(check_rising_factorial PRIVATE .)
target_link_libraries(check_rising_factorial ${EMSR_LIBRARIES} quadmath)
run_check(check_rising_factorial)

add_executable(check_shi src/check_shi.cc)
target_include_directories(check_shi PRIVATE .)
target_link_libraries(check_shi ${EMSR_LIBRARIES} quadmath)
run_check(check_shi)

add_executable(check_sinc src/check_sinc.cc)
target_include_directories(check_sinc PRIVATE .)
target_link_libraries(check_sinc ${EMSR_LIBRARIES} quadmath)
run_check(check_sinc)

add_executable(check_sinc_pi src/check_sinc_pi.cc)
target_include_directories(check_sinc_pi PRIVATE .)
target_link_libraries(check_sinc_pi ${EMSR_LIBRARIES} quadmath)
run_check(check_sinc_pi)

add_executable(check_sinhint src/check_sinhint.cc)
target_include_directories(check_sinhint PRIVATE .)
target_link_libraries(check_sinhint ${EMSR_LIBRARIES} quadmath)
run_check(check_sinhint)

add_executable(check_sinint src/check_sinint.cc)
target_include_directories(check_sinint PRIVATE .)
target_link_libraries(check_sinint ${EMSR_LIBRARIES} quadmath)
run_check(check_sinint)

add_executable(check_sin_pi src/check_sin_pi.cc)
target_include_directories(check_sin_pi PRIVATE .)
target_link_libraries(check_sin_pi ${EMSR_LIBRARIES} quadmath)
run_check(check_sin_pi)

add_executable(check_sph_bessel src/check_sph_bessel.cc)
target_include_directories(check_sph_bessel PRIVATE .)
target_link_libraries(check_sph_bessel ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_bessel)

add_executable(check_sph_bessel_i src/check_sph_bessel_i.cc)
target_include_directories(check_sph_bessel_i PRIVATE .)
target_link_libraries(check_sph_bessel_i ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_bessel_i)

add_executable(check_sph_bessel_k src/check_sph_bessel_k.cc)
target_include_directories(check_sph_bessel_k PRIVATE .)
target_link_libraries(check_sph_bessel_k ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_bessel_k)

add_executable(check_sph_hankel_1 src/check_sph_hankel_1.cc)
target_include_directories(check_sph_hankel_1 PRIVATE .)
target_link_libraries(check_sph_hankel_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_hankel_1)

add_executable(check_sph_hankel_2 src/check_sph_hankel_2.cc)
target_include_directories(check_sph_hankel_2 PRIVATE .)
target_link_libraries(check_sph_hankel_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_hankel_2)

add_executable(check_sph_harmonic src/check_sph_harmonic.cc)
target_include_directories(check_sph_harmonic PRIVATE .)
target_link_libraries(check_sph_harmonic ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_harmonic)

add_executable(check_sph_legendre src/check_sph_legendre.cc)
target_include_directories(check_sph_legendre PRIVATE .)
target_link_libraries(check_sph_legendre ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_legendre)

add_executable(check_sph_neumann src/check_sph_neumann.cc)
target_include_directories(check_sph_neumann PRIVATE .)
target_link_libraries(check_sph_neumann ${EMSR_LIBRARIES} quadmath)
run_check(check_sph_neumann)

add_executable(check_stirling_1 src/check_stirling_1.cc)
target_include_directories(check_stirling_1 PRIVATE .)
target_link_libraries(check_stirling_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_stirling_1)

add_executable(check_stirling_2 src/check_stirling_2.cc)
target_include_directories(check_stirling_2 PRIVATE .)
target_link_libraries(check_stirling_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_stirling_2)

add_executable(check_tgamma src/check_tgamma.cc)
target_include_directories(check_tgamma PRIVATE .)
target_link_libraries(check_tgamma ${EMSR_LIBRARIES} quadmath)
run_check(check_tgamma)

add_executable(check_tgamma_lower src/check_tgamma_lower.cc)
target_include_directories(check_tgamma_lower PRIVATE .)
target_link_libraries(check_tgamma_lower ${EMSR_LIBRARIES} quadmath)
run_check(check_tgamma_lower)

add_executable(check_theta_1 src/check_theta_1.cc)
target_include_directories(check_theta_1 PRIVATE .)
target_link_libraries(check_theta_1 ${EMSR_LIBRARIES} quadmath)
run_check(check_theta_1)

add_executable(check_theta_2 src/check_theta_2.cc)
target_include_directories(check_theta_2 PRIVATE .)
target_link_libraries(check_theta_2 ${EMSR_LIBRARIES} quadmath)
run_check(check_theta_2)

add_executable(check_theta_3 src/check_theta_3.cc)
target_include_directories(check_theta_3 PRIVATE .)
target_link_libraries(check_theta_3 ${EMSR_LIBRARIES} quadmath)
run_check(check_theta_3)

add_executable(check_theta_4 src/check_theta_4.cc)
target_include_directories(check_theta_4 PRIVATE .)
target_link_libraries(check_theta_4 ${EMSR_LIBRARIES} quadmath)
run_check(check_theta_4)

add_executable(check_zernike src/check_zernike.cc)
target_include_directories(check_zernike PRIVATE .)
target_link_libraries(check_zernike ${EMSR_LIBRARIES} quadmath)
run_check(check_zernike)

add_executable(complex_airy_ai src/complex_airy_ai.cc)
target_include_directories(complex_airy_ai PRIVATE .)
target_link_libraries(complex_airy_ai ${EMSR_LIBRARIES} quadmath)
run_check(complex_airy_ai)

add_executable(complex_airy_bi src/complex_airy_bi.cc)
target_include_directories(complex_airy_bi PRIVATE .)
target_link_libraries(complex_airy_bi ${EMSR_LIBRARIES} quadmath)
run_check(complex_airy_bi)

add_executable(complex_ellint_rc src/complex_ellint_rc.cc)
target_include_directories(complex_ellint_rc PRIVATE .)
target_link_libraries(complex_ellint_rc ${EMSR_LIBRARIES} quadmath)
run_check(complex_ellint_rc)

add_executable(complex_ellint_rd src/complex_ellint_rd.cc)
target_include_directories(complex_ellint_rd PRIVATE .)
target_link_libraries(complex_ellint_rd ${EMSR_LIBRARIES} quadmath)
run_check(complex_ellint_rd)

add_executable(complex_ellint_rf src/complex_ellint_rf.cc)
target_include_directories(complex_ellint_rf PRIVATE .)
target_link_libraries(complex_ellint_rf ${EMSR_LIBRARIES} quadmath)
run_check(complex_ellint_rf)

add_executable(complex_ellint_rg src/complex_ellint_rg.cc)
target_include_directories(complex_ellint_rg PRIVATE .)
target_link_libraries(complex_ellint_rg ${EMSR_LIBRARIES} quadmath)
run_check(complex_ellint_rg)

add_executable(complex_ellint_rj src/complex_ellint_rj.cc)
target_include_directories(complex_ellint_rj PRIVATE .)
target_link_libraries(complex_ellint_rj ${EMSR_LIBRARIES} quadmath)
run_check(complex_ellint_rj)

add_executable(deathmatch_comp_ellint src/deathmatch_comp_ellint.cc)
target_include_directories(deathmatch_comp_ellint PRIVATE .)
target_link_libraries(deathmatch_comp_ellint ${EMSR_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg src/deathmatch_conf_hyperg.cc)
target_include_directories(deathmatch_conf_hyperg PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg ${EMSR_LIBRARIES} quadmath)

add_executable(deathmatch_conf_hyperg_lim src/deathmatch_conf_hyperg_lim.cc)
target_include_directories(deathmatch_conf_hyperg_lim PRIVATE .)
target_link_libraries(deathmatch_conf_hyperg_lim ${EMSR_LIBRARIES} quadmath)
run_check(deathmatch_conf_hyperg_lim)

add_executable(deathmatch_hyperg src/deathmatch_hyperg.cc)
target_include_directories(deathmatch_hyperg PRIVATE .)
target_link_libraries(deathmatch_hyperg ${EMSR_LIBRARIES} quadmath)
run_check(deathmatch_hyperg)

add_executable(origin_cyl_bessel_j src/origin_cyl_bessel_j.cc)
target_include_directories(origin_cyl_bessel_j PRIVATE .)
target_link_libraries(origin_cyl_bessel_j ${EMSR_LIBRARIES} quadmath)
run_check(origin_cyl_bessel_j)

add_executable(origin_cyl_neumann src/origin_cyl_neumann.cc)
target_include_directories(origin_cyl_neumann PRIVATE .)
target_link_libraries(origin_cyl_neumann ${EMSR_LIBRARIES} quadmath)
run_check(origin_cyl_neumann)

add_executable(pr56216_cyl_bessel_i src/pr56216_cyl_bessel_i.cc)
target_include_directories(pr56216_cyl_bessel_i PRIVATE .)
target_link_libraries(pr56216_cyl_bessel_i ${EMSR_LIBRARIES} quadmath)
run_check(pr56216_cyl_bessel_i)

add_executable(pr56216_cyl_hankel_1 src/pr56216_cyl_hankel_1.cc)
target_include_directories(pr56216_cyl_hankel_1 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_1 ${EMSR_LIBRARIES} quadmath)
run_check(pr56216_cyl_hankel_1)

add_executable(pr56216_cyl_hankel_2 src/pr56216_cyl_hankel_2.cc)
target_include_directories(pr56216_cyl_hankel_2 PRIVATE .)
target_link_libraries(pr56216_cyl_hankel_2 ${EMSR_LIBRARIES} quadmath)
run_check(pr56216_cyl_hankel_2)

add_executable(pr68397 src/pr68397.cc)
target_include_directories(pr68397 PRIVATE .)
target_link_libraries(pr68397 ${EMSR_LIBRARIES} quadmath)
run_check(pr68397)

add_executable(pr83140 src/pr83140.cc)
target_include_directories(pr83140 PRIVATE .)
target_link_libraries(pr83140 ${EMSR_LIBRARIES} quadmath)
run_check(pr83140)

add_executable(pr86655_assoc_legendre src/pr86655_assoc_legendre.cc)
target_include_directories(pr86655_assoc_legendre PRIVATE .)
target_link_libraries(pr86655_assoc_legendre ${EMSR_LIBRARIES} quadmath)
run_check(pr86655_assoc_legendre)

add_executable(pr86655_sph_legendre src/pr86655_sph_legendre.cc)
target_include_directories(pr86655_sph_legendre PRIVATE .)
target_link_libraries(pr86655_sph_legendre ${EMSR_LIBRARIES} quadmath)
run_check(pr86655_sph_legendre)

# Orthonormality checks

add_executable(orthonorm_assoc_laguerre src/orthonorm_assoc_laguerre.cc)
target_include_directories(orthonorm_assoc_laguerre PRIVATE .)
target_link_libraries(orthonorm_assoc_laguerre ${EMSR_LIBRARIES})
run_check(orthonorm_assoc_laguerre)

add_executable(orthonorm_assoc_legendre src/orthonorm_assoc_legendre.cc)
target_include_directories(orthonorm_assoc_legendre PRIVATE .)
target_link_libraries(orthonorm_assoc_legendre ${EMSR_LIBRARIES})
run_check(orthonorm_assoc_legendre)

add_executable(orthonorm_chebyshev_t src/orthonorm_chebyshev_t.cc)
target_include_directories(orthonorm_chebyshev_t PRIVATE .)
target_link_libraries(orthonorm_chebyshev_t ${EMSR_LIBRARIES})
run_check(orthonorm_chebyshev_t)

add_executable(orthonorm_chebyshev_u src/orthonorm_chebyshev_u.cc)
target_include_directories(orthonorm_chebyshev_u PRIVATE .)
target_link_libraries(orthonorm_chebyshev_u ${EMSR_LIBRARIES})
run_check(orthonorm_chebyshev_u)

add_executable(orthonorm_chebyshev_v src/orthonorm_chebyshev_v.cc)
target_include_directories(orthonorm_chebyshev_v PRIVATE .)
target_link_libraries(orthonorm_chebyshev_v ${EMSR_LIBRARIES})
run_check(orthonorm_chebyshev_v)

add_executable(orthonorm_chebyshev_w src/orthonorm_chebyshev_w.cc)
target_include_directories(orthonorm_chebyshev_w PRIVATE .)
target_link_libraries(orthonorm_chebyshev_w ${EMSR_LIBRARIES})
run_check(orthonorm_chebyshev_w)

add_executable(orthonorm_gegenbauer src/orthonorm_gegenbauer.cc)
target_include_directories(orthonorm_gegenbauer PRIVATE .)
target_link_libraries(orthonorm_gegenbauer ${EMSR_LIBRARIES})
run_check(orthonorm_gegenbauer)

add_executable(orthonorm_hermite src/orthonorm_hermite.cc)
target_include_directories(orthonorm_hermite PRIVATE .)
target_link_libraries(orthonorm_hermite ${EMSR_LIBRARIES})
run_check(orthonorm_hermite)

add_executable(orthonorm_jacobi src/orthonorm_jacobi.cc)
target_include_directories(orthonorm_jacobi PRIVATE .)
target_link_libraries(orthonorm_jacobi ${EMSR_LIBRARIES})
run_check(orthonorm_jacobi)

add_executable(orthonorm_laguerre src/orthonorm_laguerre.cc)
target_include_directories(orthonorm_laguerre PRIVATE .)
target_link_libraries(orthonorm_laguerre ${EMSR_LIBRARIES})
run_check(orthonorm_laguerre)

add_executable(orthonorm_legendre src/orthonorm_legendre.cc)
target_include_directories(orthonorm_legendre PRIVATE .)
target_link_libraries(orthonorm_legendre ${EMSR_LIBRARIES})
run_check(orthonorm_legendre)

add_executable(orthonorm_radpoly src/orthonorm_radpoly.cc)
target_include_directories(orthonorm_radpoly PRIVATE .)
target_link_libraries(orthonorm_radpoly ${EMSR_LIBRARIES})
run_check(orthonorm_radpoly)

add_executable(orthonorm_sph_legendre src/orthonorm_sph_legendre.cc)
target_include_directories(orthonorm_sph_legendre PRIVATE .)
target_link_libraries(orthonorm_sph_legendre ${EMSR_LIBRARIES})
run_check(orthonorm_sph_legendre)

add_executable(orthonorm_zernike src/orthonorm_zernike.cc)
target_include_directories(orthonorm_zernike PRIVATE .)
target_link_libraries(orthonorm_zernike ${EMSR_LIBRARIES})
run_check(orthonorm_zernike)
