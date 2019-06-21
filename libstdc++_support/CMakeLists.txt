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

set(tr29124_INCLUDE_DIRS
  .
  ../wrappers
  ../include
  ../cxx_math_constants/include
  ../cxx_fp_utils/include
  ../cxx_complex_utils/include
  ../cxx_complex_math/include
  ../polynomial/include
  ../cxx_summation/include
  ../quadrature/include)

add_executable(testcase testcase.cpp)
target_link_libraries(testcase -L../wrappers quadmath -lwrap_gsl -lwrap_boost -lwrap_burkhardt)
target_include_directories(testcase PRIVATE
  ${tr29124_INCLUDE_DIRS})

add_executable(testcase_tr1 testcase.cpp)
target_link_libraries(testcase_tr1 quadmath -L../wrappers -lwrap_gsl -lwrap_boost)
target_compile_definitions(testcase_tr1 PRIVATE TR1)
target_include_directories(testcase_tr1 PRIVATE
  ${tr29124_INCLUDE_DIRS})

#
# Function to run checks.
#
function(run_check target)
  #COMMAND(bash -c "echo ${target} >> output/check_out.txt")
  COMMAND(bash -c ">&2 echo ${target} >> output/check_err.txt")
  COMMAND(bash -c "${CMAKE_BINARY_DIR}/bin/${target} >> output/check_out.txt 2>> output/check_err.txt")
endfunction()

#(target_link_libraries\(([a-z0-9_]*) quadmath\))
#\1\nadd_test(NAME run_\2 COMMAND run_check(\2))

add_executable(check_airy_ai check/check_airy_ai.cc)
target_include_directories(check_airy_ai PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_airy_ai quadmath)
add_test(NAME run_check_airy_ai COMMAND run_check(check_airy_ai))

add_executable(check_airy_bi check/check_airy_bi.cc)
target_include_directories(check_airy_bi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_airy_bi quadmath)

add_executable(check_assoc_laguerre check/check_assoc_laguerre.cc)
target_include_directories(check_assoc_laguerre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_assoc_laguerre quadmath)

add_executable(check_assoc_legendre check/check_assoc_legendre.cc)
target_include_directories(check_assoc_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_assoc_legendre quadmath)

add_executable(check_bell check/check_bell.cc)
target_include_directories(check_bell PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_bell quadmath)

add_executable(check_bernoulli check/check_bernoulli.cc)
target_include_directories(check_bernoulli PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_bernoulli quadmath)

add_executable(check_beta check/check_beta.cc)
target_include_directories(check_beta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_beta quadmath)

add_executable(check_binomial check/check_binomial.cc)
target_include_directories(check_binomial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_binomial quadmath)

add_executable(check_chebyshev_t check/check_chebyshev_t.cc)
target_include_directories(check_chebyshev_t PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_chebyshev_t quadmath)

add_executable(check_chebyshev_u check/check_chebyshev_u.cc)
target_include_directories(check_chebyshev_u PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_chebyshev_u quadmath)

add_executable(check_chebyshev_v check/check_chebyshev_v.cc)
target_include_directories(check_chebyshev_v PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_chebyshev_v quadmath)

add_executable(check_chebyshev_w check/check_chebyshev_w.cc)
target_include_directories(check_chebyshev_w PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_chebyshev_w quadmath)

add_executable(check_chi check/check_chi.cc)
target_include_directories(check_chi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_chi quadmath)

add_executable(check_clausen_cl check/check_clausen_cl.cc)
target_include_directories(check_clausen_cl PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_clausen_cl quadmath)

add_executable(check_comp_ellint_1 check/check_comp_ellint_1.cc)
target_include_directories(check_comp_ellint_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_comp_ellint_1 quadmath)

add_executable(check_comp_ellint_2 check/check_comp_ellint_2.cc)
target_include_directories(check_comp_ellint_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_comp_ellint_2 quadmath)

add_executable(check_comp_ellint_3 check/check_comp_ellint_3.cc)
target_include_directories(check_comp_ellint_3 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_comp_ellint_3 quadmath)

add_executable(check_comp_ellint_d check/check_comp_ellint_d.cc)
target_include_directories(check_comp_ellint_d PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_comp_ellint_d quadmath)

add_executable(check_conf_hyperg check/check_conf_hyperg.cc)
target_include_directories(check_conf_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_conf_hyperg quadmath)

add_executable(check_conf_hyperg_lim check/check_conf_hyperg_lim.cc)
target_include_directories(check_conf_hyperg_lim PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_conf_hyperg_lim quadmath)

add_executable(check_coshint check/check_coshint.cc)
target_include_directories(check_coshint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_coshint quadmath)

add_executable(check_cosint check/check_cosint.cc)
target_include_directories(check_cosint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cosint quadmath)

add_executable(check_cos_pi check/check_cos_pi.cc)
target_include_directories(check_cos_pi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cos_pi quadmath)

add_executable(check_cyl_bessel_i check/check_cyl_bessel_i.cc)
target_include_directories(check_cyl_bessel_i PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_bessel_i quadmath)

add_executable(check_cyl_bessel_i_scaled check/check_cyl_bessel_i_scaled.cc)
target_include_directories(check_cyl_bessel_i_scaled PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_bessel_i_scaled quadmath)

add_executable(check_cyl_bessel_j check/check_cyl_bessel_j.cc)
target_include_directories(check_cyl_bessel_j PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_bessel_j quadmath)

add_executable(check_cyl_bessel_k check/check_cyl_bessel_k.cc)
target_include_directories(check_cyl_bessel_k PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_bessel_k quadmath)

add_executable(check_cyl_bessel_k_scaled check/check_cyl_bessel_k_scaled.cc)
target_include_directories(check_cyl_bessel_k_scaled PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_bessel_k_scaled quadmath)

add_executable(check_cyl_hankel_1 check/check_cyl_hankel_1.cc)
target_include_directories(check_cyl_hankel_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_hankel_1 quadmath)

add_executable(check_cyl_hankel_2 check/check_cyl_hankel_2.cc)
target_include_directories(check_cyl_hankel_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_hankel_2 quadmath)

add_executable(check_cyl_neumann check/check_cyl_neumann.cc)
target_include_directories(check_cyl_neumann PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_cyl_neumann quadmath)

add_executable(check_dawson check/check_dawson.cc)
target_include_directories(check_dawson PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_dawson quadmath)

add_executable(check_debye check/check_debye.cc)
target_include_directories(check_debye PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_debye quadmath)

add_executable(check_digamma check/check_digamma.cc)
target_include_directories(check_digamma PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_digamma quadmath)

add_executable(check_dilog check/check_dilog.cc)
target_include_directories(check_dilog PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_dilog quadmath)

add_executable(check_dirichlet_beta check/check_dirichlet_beta.cc)
target_include_directories(check_dirichlet_beta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_dirichlet_beta quadmath)

add_executable(check_dirichlet_eta check/check_dirichlet_eta.cc)
target_include_directories(check_dirichlet_eta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_dirichlet_eta quadmath)

add_executable(check_dirichlet_lambda check/check_dirichlet_lambda.cc)
target_include_directories(check_dirichlet_lambda PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_dirichlet_lambda quadmath)

add_executable(check_double_factorial check/check_double_factorial.cc)
target_include_directories(check_double_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_double_factorial quadmath)

add_executable(check_ellint_1 check/check_ellint_1.cc)
target_include_directories(check_ellint_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_1 quadmath)

add_executable(check_ellint_2 check/check_ellint_2.cc)
target_include_directories(check_ellint_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_2 quadmath)

add_executable(check_ellint_3 check/check_ellint_3.cc)
target_include_directories(check_ellint_3 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_3 quadmath)

add_executable(check_ellint_d check/check_ellint_d.cc)
target_include_directories(check_ellint_d PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_d quadmath)

add_executable(check_ellint_rc check/check_ellint_rc.cc)
target_include_directories(check_ellint_rc PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_rc quadmath)

add_executable(check_ellint_rd check/check_ellint_rd.cc)
target_include_directories(check_ellint_rd PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_rd quadmath)

add_executable(check_ellint_rf check/check_ellint_rf.cc)
target_include_directories(check_ellint_rf PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_rf quadmath)

add_executable(check_ellint_rg check/check_ellint_rg.cc)
target_include_directories(check_ellint_rg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_rg quadmath)

add_executable(check_ellint_rj check/check_ellint_rj.cc)
target_include_directories(check_ellint_rj PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellint_rj quadmath)

add_executable(check_ellnome check/check_ellnome.cc)
target_include_directories(check_ellnome PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ellnome quadmath)

add_executable(check_euler check/check_euler.cc)
target_include_directories(check_euler PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_euler quadmath)

add_executable(check_eulerian_1 check/check_eulerian_1.cc)
target_include_directories(check_eulerian_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_eulerian_1 quadmath)

add_executable(check_eulerian_2 check/check_eulerian_2.cc)
target_include_directories(check_eulerian_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_eulerian_2 quadmath)

add_executable(check_expint check/check_expint.cc)
target_include_directories(check_expint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_expint quadmath)

add_executable(check_expint_en check/check_expint_en.cc)
target_include_directories(check_expint_en PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_expint_en quadmath)

add_executable(check_factorial check/check_factorial.cc)
target_include_directories(check_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_factorial quadmath)

add_executable(check_falling_factorial check/check_falling_factorial.cc)
target_include_directories(check_falling_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_falling_factorial quadmath)

add_executable(check_fresnel_c check/check_fresnel_c.cc)
target_include_directories(check_fresnel_c PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_fresnel_c quadmath)

add_executable(check_fresnel_s check/check_fresnel_s.cc)
target_include_directories(check_fresnel_s PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_fresnel_s quadmath)

add_executable(check_gamma_p check/check_gamma_p.cc)
target_include_directories(check_gamma_p PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_gamma_p quadmath)

add_executable(check_gamma_q check/check_gamma_q.cc)
target_include_directories(check_gamma_q PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_gamma_q quadmath)

add_executable(check_gamma_reciprocal check/check_gamma_reciprocal.cc)
target_include_directories(check_gamma_reciprocal PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_gamma_reciprocal quadmath)

add_executable(check_gegenbauer check/check_gegenbauer.cc)
target_include_directories(check_gegenbauer PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_gegenbauer quadmath)

add_executable(check_hermite check/check_hermite.cc)
target_include_directories(check_hermite PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_hermite quadmath)

add_executable(check_heuman_lambda check/check_heuman_lambda.cc)
target_include_directories(check_heuman_lambda PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_heuman_lambda quadmath)

add_executable(check_hurwitz_zeta check/check_hurwitz_zeta.cc)
target_include_directories(check_hurwitz_zeta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_hurwitz_zeta quadmath)

add_executable(check_hyperg check/check_hyperg.cc)
target_include_directories(check_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_hyperg quadmath)

add_executable(check_ibeta check/check_ibeta.cc)
target_include_directories(check_ibeta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ibeta quadmath)

add_executable(check_ibetac check/check_ibetac.cc)
target_include_directories(check_ibetac PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ibetac quadmath)

add_executable(check_jacobi check/check_jacobi.cc)
target_include_directories(check_jacobi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_jacobi quadmath)

add_executable(check_jacobi_cn check/check_jacobi_cn.cc)
target_include_directories(check_jacobi_cn PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_jacobi_cn quadmath)

add_executable(check_jacobi_dn check/check_jacobi_dn.cc)
target_include_directories(check_jacobi_dn PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_jacobi_dn quadmath)

add_executable(check_jacobi_sn check/check_jacobi_sn.cc)
target_include_directories(check_jacobi_sn PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_jacobi_sn quadmath)

add_executable(check_jacobi_zeta check/check_jacobi_zeta.cc)
target_include_directories(check_jacobi_zeta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_jacobi_zeta quadmath)

add_executable(check_laguerre check/check_laguerre.cc)
target_include_directories(check_laguerre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_laguerre quadmath)

add_executable(check_lah check/check_lah.cc)
target_include_directories(check_lah PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lah quadmath)

add_executable(check_lbinomial check/check_lbinomial.cc)
target_include_directories(check_lbinomial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lbinomial quadmath)

add_executable(check_ldouble_factorial check/check_ldouble_factorial.cc)
target_include_directories(check_ldouble_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_ldouble_factorial quadmath)

add_executable(check_legendre check/check_legendre.cc)
target_include_directories(check_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_legendre quadmath)

add_executable(check_legendre_q check/check_legendre_q.cc)
target_include_directories(check_legendre_q PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_legendre_q quadmath)

add_executable(check_lfactorial check/check_lfactorial.cc)
target_include_directories(check_lfactorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lfactorial quadmath)

add_executable(check_lfalling_factorial check/check_lfalling_factorial.cc)
target_include_directories(check_lfalling_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lfalling_factorial quadmath)

add_executable(check_lgamma check/check_lgamma.cc)
target_include_directories(check_lgamma PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lgamma quadmath)

add_executable(check_logistic_p check/check_logistic_p.cc)
target_include_directories(check_logistic_p PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_logistic_p quadmath)

add_executable(check_logistic_pdf check/check_logistic_pdf.cc)
target_include_directories(check_logistic_pdf PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_logistic_pdf quadmath)

add_executable(check_lognormal_p check/check_lognormal_p.cc)
target_include_directories(check_lognormal_p PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lognormal_p quadmath)

add_executable(check_lognormal_pdf check/check_lognormal_pdf.cc)
target_include_directories(check_lognormal_pdf PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lognormal_pdf quadmath)

add_executable(check_lrising_factorial check/check_lrising_factorial.cc)
target_include_directories(check_lrising_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_lrising_factorial quadmath)

add_executable(check_normal_p check/check_normal_p.cc)
target_include_directories(check_normal_p PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_normal_p quadmath)

add_executable(check_normal_pdf check/check_normal_pdf.cc)
target_include_directories(check_normal_pdf PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_normal_pdf quadmath)

add_executable(check_owens_t check/check_owens_t.cc)
target_include_directories(check_owens_t PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_owens_t quadmath)

add_executable(check_polygamma check/check_polygamma.cc)
target_include_directories(check_polygamma PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_polygamma quadmath)

add_executable(check_radpoly check/check_radpoly.cc)
target_include_directories(check_radpoly PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_radpoly quadmath)

add_executable(check_riemann_zeta check/check_riemann_zeta.cc)
target_include_directories(check_riemann_zeta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_riemann_zeta quadmath)

add_executable(check_rising_factorial check/check_rising_factorial.cc)
target_include_directories(check_rising_factorial PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_rising_factorial quadmath)

add_executable(check_shi check/check_shi.cc)
target_include_directories(check_shi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_shi quadmath)

add_executable(check_sinc check/check_sinc.cc)
target_include_directories(check_sinc PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sinc quadmath)

add_executable(check_sinc_pi check/check_sinc_pi.cc)
target_include_directories(check_sinc_pi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sinc_pi quadmath)

add_executable(check_sinhint check/check_sinhint.cc)
target_include_directories(check_sinhint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sinhint quadmath)

add_executable(check_sinint check/check_sinint.cc)
target_include_directories(check_sinint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sinint quadmath)

add_executable(check_sin_pi check/check_sin_pi.cc)
target_include_directories(check_sin_pi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sin_pi quadmath)

add_executable(check_sph_bessel check/check_sph_bessel.cc)
target_include_directories(check_sph_bessel PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_bessel quadmath)

add_executable(check_sph_bessel_i check/check_sph_bessel_i.cc)
target_include_directories(check_sph_bessel_i PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_bessel_i quadmath)

add_executable(check_sph_bessel_k check/check_sph_bessel_k.cc)
target_include_directories(check_sph_bessel_k PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_bessel_k quadmath)

add_executable(check_sph_hankel_1 check/check_sph_hankel_1.cc)
target_include_directories(check_sph_hankel_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_hankel_1 quadmath)

add_executable(check_sph_hankel_2 check/check_sph_hankel_2.cc)
target_include_directories(check_sph_hankel_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_hankel_2 quadmath)

add_executable(check_sph_harmonic check/check_sph_harmonic.cc)
target_include_directories(check_sph_harmonic PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_harmonic quadmath)

add_executable(check_sph_legendre check/check_sph_legendre.cc)
target_include_directories(check_sph_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_legendre quadmath)

add_executable(check_sph_neumann check/check_sph_neumann.cc)
target_include_directories(check_sph_neumann PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_sph_neumann quadmath)

add_executable(check_stirling_1 check/check_stirling_1.cc)
target_include_directories(check_stirling_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_stirling_1 quadmath)

add_executable(check_stirling_2 check/check_stirling_2.cc)
target_include_directories(check_stirling_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_stirling_2 quadmath)

add_executable(check_tgamma check/check_tgamma.cc)
target_include_directories(check_tgamma PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tgamma quadmath)

add_executable(check_tgamma_lower check/check_tgamma_lower.cc)
target_include_directories(check_tgamma_lower PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tgamma_lower quadmath)

add_executable(check_theta_1 check/check_theta_1.cc)
target_include_directories(check_theta_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_theta_1 quadmath)

add_executable(check_theta_2 check/check_theta_2.cc)
target_include_directories(check_theta_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_theta_2 quadmath)

add_executable(check_theta_3 check/check_theta_3.cc)
target_include_directories(check_theta_3 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_theta_3 quadmath)

add_executable(check_theta_4 check/check_theta_4.cc)
target_include_directories(check_theta_4 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_theta_4 quadmath)

add_executable(check_tr1_assoc_laguerre check/check_tr1_assoc_laguerre.cc)
target_include_directories(check_tr1_assoc_laguerre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_assoc_laguerre quadmath)

add_executable(check_tr1_assoc_legendre check/check_tr1_assoc_legendre.cc)
target_include_directories(check_tr1_assoc_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_assoc_legendre quadmath)

add_executable(check_tr1_beta check/check_tr1_beta.cc)
target_include_directories(check_tr1_beta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_beta quadmath)

add_executable(check_tr1_comp_ellint_1 check/check_tr1_comp_ellint_1.cc)
target_include_directories(check_tr1_comp_ellint_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_comp_ellint_1 quadmath)

add_executable(check_tr1_comp_ellint_2 check/check_tr1_comp_ellint_2.cc)
target_include_directories(check_tr1_comp_ellint_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_comp_ellint_2 quadmath)

add_executable(check_tr1_comp_ellint_3 check/check_tr1_comp_ellint_3.cc)
target_include_directories(check_tr1_comp_ellint_3 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_comp_ellint_3 quadmath)

add_executable(check_tr1_conf_hyperg check/check_tr1_conf_hyperg.cc)
target_include_directories(check_tr1_conf_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_conf_hyperg quadmath)

add_executable(check_tr1_cyl_bessel_i check/check_tr1_cyl_bessel_i.cc)
target_include_directories(check_tr1_cyl_bessel_i PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_cyl_bessel_i quadmath)

add_executable(check_tr1_cyl_bessel_j check/check_tr1_cyl_bessel_j.cc)
target_include_directories(check_tr1_cyl_bessel_j PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_cyl_bessel_j quadmath)

add_executable(check_tr1_cyl_bessel_k check/check_tr1_cyl_bessel_k.cc)
target_include_directories(check_tr1_cyl_bessel_k PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_cyl_bessel_k quadmath)

add_executable(check_tr1_cyl_neumann check/check_tr1_cyl_neumann.cc)
target_include_directories(check_tr1_cyl_neumann PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_cyl_neumann quadmath)

add_executable(check_tr1_ellint_1 check/check_tr1_ellint_1.cc)
target_include_directories(check_tr1_ellint_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_ellint_1 quadmath)

add_executable(check_tr1_ellint_2 check/check_tr1_ellint_2.cc)
target_include_directories(check_tr1_ellint_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_ellint_2 quadmath)

add_executable(check_tr1_ellint_3 check/check_tr1_ellint_3.cc)
target_include_directories(check_tr1_ellint_3 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_ellint_3 quadmath)

add_executable(check_tr1_expint check/check_tr1_expint.cc)
target_include_directories(check_tr1_expint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_expint quadmath)

add_executable(check_tr1_hermite check/check_tr1_hermite.cc)
target_include_directories(check_tr1_hermite PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_hermite quadmath)

add_executable(check_tr1_hyperg check/check_tr1_hyperg.cc)
target_include_directories(check_tr1_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_hyperg quadmath)

add_executable(check_tr1_laguerre check/check_tr1_laguerre.cc)
target_include_directories(check_tr1_laguerre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_laguerre quadmath)

add_executable(check_tr1_legendre check/check_tr1_legendre.cc)
target_include_directories(check_tr1_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_legendre quadmath)

add_executable(check_tr1_riemann_zeta check/check_tr1_riemann_zeta.cc)
target_include_directories(check_tr1_riemann_zeta PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_riemann_zeta quadmath)

add_executable(check_tr1_sph_bessel check/check_tr1_sph_bessel.cc)
target_include_directories(check_tr1_sph_bessel PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_sph_bessel quadmath)

add_executable(check_tr1_sph_legendre check/check_tr1_sph_legendre.cc)
target_include_directories(check_tr1_sph_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_sph_legendre quadmath)

add_executable(check_tr1_sph_neumann check/check_tr1_sph_neumann.cc)
target_include_directories(check_tr1_sph_neumann PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_tr1_sph_neumann quadmath)

add_executable(check_zernike check/check_zernike.cc)
target_include_directories(check_zernike PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(check_zernike quadmath)

add_executable(complex_airy_ai check/complex_airy_ai.cc)
target_include_directories(complex_airy_ai PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_airy_ai quadmath)

add_executable(complex_airy_bi check/complex_airy_bi.cc)
target_include_directories(complex_airy_bi PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_airy_bi quadmath)

add_executable(complex_ellint_rc check/complex_ellint_rc.cc)
target_include_directories(complex_ellint_rc PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_ellint_rc quadmath)

add_executable(complex_ellint_rd check/complex_ellint_rd.cc)
target_include_directories(complex_ellint_rd PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_ellint_rd quadmath)

add_executable(complex_ellint_rf check/complex_ellint_rf.cc)
target_include_directories(complex_ellint_rf PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_ellint_rf quadmath)

add_executable(complex_ellint_rg check/complex_ellint_rg.cc)
target_include_directories(complex_ellint_rg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_ellint_rg quadmath)

add_executable(complex_ellint_rj check/complex_ellint_rj.cc)
target_include_directories(complex_ellint_rj PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(complex_ellint_rj quadmath)

add_executable(deathmatch_comp_ellint check/deathmatch_comp_ellint.cc)
target_include_directories(deathmatch_comp_ellint PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(deathmatch_comp_ellint quadmath)

add_executable(deathmatch_conf_hyperg check/deathmatch_conf_hyperg.cc)
target_include_directories(deathmatch_conf_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(deathmatch_conf_hyperg quadmath)

add_executable(deathmatch_conf_hyperg_lim check/deathmatch_conf_hyperg_lim.cc)
target_include_directories(deathmatch_conf_hyperg_lim PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(deathmatch_conf_hyperg_lim quadmath)

add_executable(deathmatch_hyperg check/deathmatch_hyperg.cc)
target_include_directories(deathmatch_hyperg PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(deathmatch_hyperg quadmath)

add_executable(origin_cyl_bessel_j check/origin_cyl_bessel_j.cc)
target_include_directories(origin_cyl_bessel_j PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(origin_cyl_bessel_j quadmath)

add_executable(origin_cyl_neumann check/origin_cyl_neumann.cc)
target_include_directories(origin_cyl_neumann PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(origin_cyl_neumann quadmath)

add_executable(pr56216_cyl_bessel_i check/pr56216_cyl_bessel_i.cc)
target_include_directories(pr56216_cyl_bessel_i PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr56216_cyl_bessel_i quadmath)

add_executable(pr56216_cyl_hankel_1 check/pr56216_cyl_hankel_1.cc)
target_include_directories(pr56216_cyl_hankel_1 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr56216_cyl_hankel_1 quadmath)

add_executable(pr56216_cyl_hankel_2 check/pr56216_cyl_hankel_2.cc)
target_include_directories(pr56216_cyl_hankel_2 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr56216_cyl_hankel_2 quadmath)

add_executable(pr68397 check/pr68397.cc)
target_include_directories(pr68397 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr68397 quadmath)

add_executable(pr83140 check/pr83140.cc)
target_include_directories(pr83140 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr83140 quadmath)

add_executable(pr86655_assoc_legendre check/pr86655_assoc_legendre.cc)
target_include_directories(pr86655_assoc_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr86655_assoc_legendre quadmath)

add_executable(pr86655_sph_legendre check/pr86655_sph_legendre.cc)
target_include_directories(pr86655_sph_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr86655_sph_legendre quadmath)

add_executable(pr86655_tr1_assoc_legendre check/pr86655_tr1_assoc_legendre.cc)
target_include_directories(pr86655_tr1_assoc_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr86655_tr1_assoc_legendre quadmath)

add_executable(pr86655_tr1_sph_legendre check/pr86655_tr1_sph_legendre.cc)
target_include_directories(pr86655_tr1_sph_legendre PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(pr86655_tr1_sph_legendre quadmath)

# This needs C++14!
#add_executable(tr1_airy check/tr1_airy.cc)
#target_include_directories(tr1_airy PRIVATE ${tr29124_INCLUDE_DIRS})
#target_link_libraries(tr1_airy quadmath)

add_executable(tr1_pr68397 check/tr1_pr68397.cc)
target_include_directories(tr1_pr68397 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(tr1_pr68397 quadmath)

add_executable(tr1_pr83140 check/tr1_pr83140.cc)
target_include_directories(tr1_pr83140 PRIVATE ${tr29124_INCLUDE_DIRS})
target_link_libraries(tr1_pr83140 quadmath)

