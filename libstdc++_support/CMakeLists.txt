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

set(TR29124_INCLUDE_DIRS
  ../include
  ../cxx_math_constants/include
  ../cxx_fp_utils/include
  ../cxx_complex_utils/include
  ../cxx_complex_math/include
  ../cxx_continued_fractions/include
  ../cxx_summation/include
  ../cxx_polynomial/include
  ../cxx_integration/include)

add_executable(testcase testcase.cpp)
target_link_libraries(testcase -L../wrappers quadmath -lwrap_gsl -lwrap_boost -lwrap_burkhardt)
target_include_directories(testcase PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)

add_executable(testcase_tr1 testcase.cpp)
target_link_libraries(testcase_tr1 quadmath -L../wrappers -lwrap_gsl -lwrap_boost)
target_compile_definitions(testcase_tr1 PRIVATE TR1)
target_include_directories(testcase_tr1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)

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
target_include_directories(check_airy_ai PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_airy_ai quadmath)
run_check(check_airy_ai)

add_executable(check_airy_bi check/check_airy_bi.cc)
target_include_directories(check_airy_bi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_airy_bi quadmath)
run_check(check_airy_bi)

add_executable(check_assoc_laguerre check/check_assoc_laguerre.cc)
target_include_directories(check_assoc_laguerre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_assoc_laguerre quadmath)
run_check(check_assoc_laguerre)

add_executable(check_assoc_legendre check/check_assoc_legendre.cc)
target_include_directories(check_assoc_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_assoc_legendre quadmath)
run_check(check_assoc_legendre)

add_executable(check_bell check/check_bell.cc)
target_include_directories(check_bell PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_bell quadmath)
run_check(check_bell)

add_executable(check_bernoulli check/check_bernoulli.cc)
target_include_directories(check_bernoulli PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_bernoulli quadmath)
run_check(check_bernoulli)

add_executable(check_beta check/check_beta.cc)
target_include_directories(check_beta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_beta quadmath)
run_check(check_beta)

add_executable(check_binomial check/check_binomial.cc)
target_include_directories(check_binomial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_binomial quadmath)
run_check(check_binomial)

add_executable(check_chebyshev_t check/check_chebyshev_t.cc)
target_include_directories(check_chebyshev_t PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_chebyshev_t quadmath)
run_check(check_chebyshev_t)

add_executable(check_chebyshev_u check/check_chebyshev_u.cc)
target_include_directories(check_chebyshev_u PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_chebyshev_u quadmath)
run_check(check_chebyshev_u)

add_executable(check_chebyshev_v check/check_chebyshev_v.cc)
target_include_directories(check_chebyshev_v PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_chebyshev_v quadmath)
run_check(check_chebyshev_v)

add_executable(check_chebyshev_w check/check_chebyshev_w.cc)
target_include_directories(check_chebyshev_w PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_chebyshev_w quadmath)
run_check(check_chebyshev_w)

add_executable(check_chi check/check_chi.cc)
target_include_directories(check_chi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_chi quadmath)
run_check(check_chi)

add_executable(check_clausen_cl check/check_clausen_cl.cc)
target_include_directories(check_clausen_cl PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_clausen_cl quadmath)
run_check(check_clausen_cl)

add_executable(check_comp_ellint_1 check/check_comp_ellint_1.cc)
target_include_directories(check_comp_ellint_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_comp_ellint_1 quadmath)
run_check(check_comp_ellint_1)

add_executable(check_comp_ellint_2 check/check_comp_ellint_2.cc)
target_include_directories(check_comp_ellint_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_comp_ellint_2 quadmath)
run_check(check_comp_ellint_2)

add_executable(check_comp_ellint_3 check/check_comp_ellint_3.cc)
target_include_directories(check_comp_ellint_3 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_comp_ellint_3 quadmath)
run_check(check_comp_ellint_3)

add_executable(check_comp_ellint_d check/check_comp_ellint_d.cc)
target_include_directories(check_comp_ellint_d PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_comp_ellint_d quadmath)
run_check(check_comp_ellint_d)

add_executable(check_conf_hyperg check/check_conf_hyperg.cc)
target_include_directories(check_conf_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_conf_hyperg quadmath)
run_check(check_conf_hyperg)

add_executable(check_conf_hyperg_lim check/check_conf_hyperg_lim.cc)
target_include_directories(check_conf_hyperg_lim PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_conf_hyperg_lim quadmath)
run_check(check_conf_hyperg_lim)

add_executable(check_coshint check/check_coshint.cc)
target_include_directories(check_coshint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_coshint quadmath)
run_check(check_coshint)

add_executable(check_cosint check/check_cosint.cc)
target_include_directories(check_cosint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cosint quadmath)
run_check(check_cosint)

add_executable(check_cos_pi check/check_cos_pi.cc)
target_include_directories(check_cos_pi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cos_pi quadmath)
run_check(check_cos_pi)

add_executable(check_cyl_bessel_i check/check_cyl_bessel_i.cc)
target_include_directories(check_cyl_bessel_i PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_bessel_i quadmath)
run_check(check_cyl_bessel_i)

add_executable(check_cyl_bessel_i_scaled check/check_cyl_bessel_i_scaled.cc)
target_include_directories(check_cyl_bessel_i_scaled PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_bessel_i_scaled quadmath)
run_check(check_cyl_bessel_i_scaled)

add_executable(check_cyl_bessel_j check/check_cyl_bessel_j.cc)
target_include_directories(check_cyl_bessel_j PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_bessel_j quadmath)
run_check(check_cyl_bessel_j)

add_executable(check_cyl_bessel_k check/check_cyl_bessel_k.cc)
target_include_directories(check_cyl_bessel_k PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_bessel_k quadmath)
run_check(check_cyl_bessel_k)

add_executable(check_cyl_bessel_k_scaled check/check_cyl_bessel_k_scaled.cc)
target_include_directories(check_cyl_bessel_k_scaled PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_bessel_k_scaled quadmath)
run_check(check_cyl_bessel_k_scaled)

add_executable(check_cyl_hankel_1 check/check_cyl_hankel_1.cc)
target_include_directories(check_cyl_hankel_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_hankel_1 quadmath)
run_check(check_cyl_hankel_1)

add_executable(check_cyl_hankel_2 check/check_cyl_hankel_2.cc)
target_include_directories(check_cyl_hankel_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_hankel_2 quadmath)
run_check(check_cyl_hankel_2)

add_executable(check_cyl_neumann check/check_cyl_neumann.cc)
target_include_directories(check_cyl_neumann PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_cyl_neumann quadmath)
run_check(check_cyl_neumann)

add_executable(check_dawson check/check_dawson.cc)
target_include_directories(check_dawson PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_dawson quadmath)
run_check(check_dawson)

add_executable(check_debye check/check_debye.cc)
target_include_directories(check_debye PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_debye quadmath)
run_check(check_debye)

add_executable(check_digamma check/check_digamma.cc)
target_include_directories(check_digamma PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_digamma quadmath)
run_check(check_digamma)

add_executable(check_dilog check/check_dilog.cc)
target_include_directories(check_dilog PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_dilog quadmath)
run_check(check_dilog)

add_executable(check_dirichlet_beta check/check_dirichlet_beta.cc)
target_include_directories(check_dirichlet_beta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_dirichlet_beta quadmath)
run_check(check_dirichlet_beta)

add_executable(check_dirichlet_eta check/check_dirichlet_eta.cc)
target_include_directories(check_dirichlet_eta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_dirichlet_eta quadmath)
run_check(check_dirichlet_eta)

add_executable(check_dirichlet_lambda check/check_dirichlet_lambda.cc)
target_include_directories(check_dirichlet_lambda PRIVATE
  .run_check(check_dirichlet_lambda)

  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_dirichlet_lambda quadmath)
run_check(check_dirichlet_lambda)

add_executable(check_double_factorial check/check_double_factorial.cc)
target_include_directories(check_double_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_double_factorial quadmath)
run_check(check_double_factorial)

add_executable(check_ellint_1 check/check_ellint_1.cc)
target_include_directories(check_ellint_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_1 quadmath)
run_check(check_ellint_1)

add_executable(check_ellint_2 check/check_ellint_2.cc)
target_include_directories(check_ellint_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_2 quadmath)
run_check(check_ellint_2)

add_executable(check_ellint_3 check/check_ellint_3.cc)
target_include_directories(check_ellint_3 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_3 quadmath)
run_check(check_ellint_3)

add_executable(check_ellint_d check/check_ellint_d.cc)
target_include_directories(check_ellint_d PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_d quadmath)

add_executable(check_ellint_rc check/check_ellint_rc.cc)
target_include_directories(check_ellint_rc PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_rc quadmath)

add_executable(check_ellint_rd check/check_ellint_rd.cc)
target_include_directories(check_ellint_rd PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_rd quadmath)
run_check(check_ellint_rd)

add_executable(check_ellint_rf check/check_ellint_rf.cc)
target_include_directories(check_ellint_rf PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_rf quadmath)
run_check(check_ellint_rf)

add_executable(check_ellint_rg check/check_ellint_rg.cc)
target_include_directories(check_ellint_rg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_rg quadmath)
run_check(check_ellint_rg)

add_executable(check_ellint_rj check/check_ellint_rj.cc)
target_include_directories(check_ellint_rj PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellint_rj quadmath)
run_check(check_ellint_rj)

add_executable(check_ellnome check/check_ellnome.cc)
target_include_directories(check_ellnome PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ellnome quadmath)
run_check(check_ellnome)

add_executable(check_euler check/check_euler.cc)
target_include_directories(check_euler PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_euler quadmath)
run_check(check_euler)

add_executable(check_eulerian_1 check/check_eulerian_1.cc)
target_include_directories(check_eulerian_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_eulerian_1 quadmath)
run_check(check_eulerian_1)

add_executable(check_eulerian_2 check/check_eulerian_2.cc)
target_include_directories(check_eulerian_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_eulerian_2 quadmath)
run_check(check_eulerian_2)

add_executable(check_expint check/check_expint.cc)
target_include_directories(check_expint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_expint quadmath)
run_check(check_expint)

add_executable(check_expint_en check/check_expint_en.cc)
target_include_directories(check_expint_en PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_expint_en quadmath)
run_check(check_expint_en)

add_executable(check_factorial check/check_factorial.cc)
target_include_directories(check_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_factorial quadmath)
run_check(check_factorial)

add_executable(check_falling_factorial check/check_falling_factorial.cc)
target_include_directories(check_falling_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_falling_factorial quadmath)
run_check(check_falling_factorial)

add_executable(check_fresnel_c check/check_fresnel_c.cc)
target_include_directories(check_fresnel_c PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_fresnel_c quadmath)
run_check(check_fresnel_c)

add_executable(check_fresnel_s check/check_fresnel_s.cc)
target_include_directories(check_fresnel_s PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_fresnel_s quadmath)
run_check(check_fresnel_s)

add_executable(check_gamma_p check/check_gamma_p.cc)
target_include_directories(check_gamma_p PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_gamma_p quadmath)
run_check(check_gamma_p)

add_executable(check_gamma_q check/check_gamma_q.cc)
target_include_directories(check_gamma_q PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_gamma_q quadmath)
run_check(check_gamma_q)

add_executable(check_gamma_reciprocal check/check_gamma_reciprocal.cc)
target_include_directories(check_gamma_reciprocal PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_gamma_reciprocal quadmath)
run_check(check_gamma_reciprocal)

add_executable(check_gegenbauer check/check_gegenbauer.cc)
target_include_directories(check_gegenbauer PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_gegenbauer quadmath)
run_check(check_gegenbauer)

add_executable(check_hermite check/check_hermite.cc)
target_include_directories(check_hermite PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_hermite quadmath)
run_check(check_hermite)

add_executable(check_heuman_lambda check/check_heuman_lambda.cc)
target_include_directories(check_heuman_lambda PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_heuman_lambda quadmath)
run_check(check_heuman_lambda)

add_executable(check_hurwitz_zeta check/check_hurwitz_zeta.cc)
target_include_directories(check_hurwitz_zeta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_hurwitz_zeta quadmath)
run_check(check_hurwitz_zeta)

add_executable(check_hyperg check/check_hyperg.cc)
target_include_directories(check_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_hyperg quadmath)
run_check(check_hyperg)

add_executable(check_ibeta check/check_ibeta.cc)
target_include_directories(check_ibeta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ibeta quadmath)
run_check(check_ibeta)

add_executable(check_ibetac check/check_ibetac.cc)
target_include_directories(check_ibetac PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ibetac quadmath)
run_check(check_ibetac)

add_executable(check_jacobi check/check_jacobi.cc)
target_include_directories(check_jacobi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_jacobi quadmath)

add_executable(check_jacobi_cn check/check_jacobi_cn.cc)
target_include_directories(check_jacobi_cn PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_jacobi_cn quadmath)
run_check(check_jacobi_cn)

add_executable(check_jacobi_dn check/check_jacobi_dn.cc)
target_include_directories(check_jacobi_dn PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_jacobi_dn quadmath)
run_check(check_jacobi_dn)

add_executable(check_jacobi_sn check/check_jacobi_sn.cc)
target_include_directories(check_jacobi_sn PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_jacobi_sn quadmath)
run_check(check_jacobi_sn)

add_executable(check_jacobi_zeta check/check_jacobi_zeta.cc)
target_include_directories(check_jacobi_zeta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_jacobi_zeta quadmath)
run_check(check_jacobi_zeta)

add_executable(check_laguerre check/check_laguerre.cc)
target_include_directories(check_laguerre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_laguerre quadmath)
run_check(check_laguerre)

add_executable(check_lah check/check_lah.cc)
target_include_directories(check_lah PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lah quadmath)
run_check(check_lah)

add_executable(check_lbinomial check/check_lbinomial.cc)
target_include_directories(check_lbinomial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lbinomial quadmath)
run_check(check_lbinomial)

add_executable(check_ldouble_factorial check/check_ldouble_factorial.cc)
target_include_directories(check_ldouble_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_ldouble_factorial quadmath)
run_check(check_ldouble_factorial)

add_executable(check_legendre check/check_legendre.cc)
target_include_directories(check_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_legendre quadmath)
run_check(check_legendre)

add_executable(check_legendre_q check/check_legendre_q.cc)
target_include_directories(check_legendre_q PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_legendre_q quadmath)
run_check(check_legendre_q)

add_executable(check_lfactorial check/check_lfactorial.cc)
target_include_directories(check_lfactorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lfactorial quadmath)
run_check(check_lfactorial)

add_executable(check_lfalling_factorial check/check_lfalling_factorial.cc)
target_include_directories(check_lfalling_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lfalling_factorial quadmath)
run_check(check_lfalling_factorial)

add_executable(check_lgamma check/check_lgamma.cc)
target_include_directories(check_lgamma PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lgamma quadmath)
run_check(check_lgamma)

add_executable(check_logistic_p check/check_logistic_p.cc)
target_include_directories(check_logistic_p PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_logistic_p quadmath)
run_check(check_logistic_p)

add_executable(check_logistic_pdf check/check_logistic_pdf.cc)
target_include_directories(check_logistic_pdf PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_logistic_pdf quadmath)
run_check(check_logistic_pdf)

add_executable(check_lognormal_p check/check_lognormal_p.cc)
target_include_directories(check_lognormal_p PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lognormal_p quadmath)
run_check(check_lognormal_p)

add_executable(check_lognormal_pdf check/check_lognormal_pdf.cc)
target_include_directories(check_lognormal_pdf PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lognormal_pdf quadmath)
run_check(check_lognormal_pdf)

add_executable(check_lrising_factorial check/check_lrising_factorial.cc)
target_include_directories(check_lrising_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_lrising_factorial quadmath)
run_check(check_lrising_factorial)

add_executable(check_normal_p check/check_normal_p.cc)
target_include_directories(check_normal_p PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_normal_p quadmath)
run_check(check_normal_p)

add_executable(check_normal_pdf check/check_normal_pdf.cc)
target_include_directories(check_normal_pdf PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_normal_pdf quadmath)
run_check(check_normal_pdf)

add_executable(check_owens_t check/check_owens_t.cc)
target_include_directories(check_owens_t PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_owens_t quadmath)
run_check(check_owens_t)

add_executable(check_polygamma check/check_polygamma.cc)
target_include_directories(check_polygamma PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_polygamma quadmath)
run_check(check_polygamma)

add_executable(check_radpoly check/check_radpoly.cc)
target_include_directories(check_radpoly PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_radpoly quadmath)
run_check(check_radpoly)

add_executable(check_riemann_zeta check/check_riemann_zeta.cc)
target_include_directories(check_riemann_zeta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_riemann_zeta quadmath)
run_check(check_riemann_zeta)

add_executable(check_rising_factorial check/check_rising_factorial.cc)
target_include_directories(check_rising_factorial PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_rising_factorial quadmath)
run_check(check_rising_factorial)

add_executable(check_shi check/check_shi.cc)
target_include_directories(check_shi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_shi quadmath)
run_check(check_shi)

add_executable(check_sinc check/check_sinc.cc)
target_include_directories(check_sinc PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sinc quadmath)
run_check(check_sinc)

add_executable(check_sinc_pi check/check_sinc_pi.cc)
target_include_directories(check_sinc_pi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sinc_pi quadmath)
run_check(check_sinc_pi)

add_executable(check_sinhint check/check_sinhint.cc)
target_include_directories(check_sinhint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sinhint quadmath)
run_check(check_sinhint)

add_executable(check_sinint check/check_sinint.cc)
target_include_directories(check_sinint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sinint quadmath)
run_check(check_sinint)

add_executable(check_sin_pi check/check_sin_pi.cc)
target_include_directories(check_sin_pi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sin_pi quadmath)
run_check(check_sin_pi)

add_executable(check_sph_bessel check/check_sph_bessel.cc)
target_include_directories(check_sph_bessel PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_bessel quadmath)
run_check(check_sph_bessel)

add_executable(check_sph_bessel_i check/check_sph_bessel_i.cc)
target_include_directories(check_sph_bessel_i PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_bessel_i quadmath)
run_check(check_sph_bessel_i)

add_executable(check_sph_bessel_k check/check_sph_bessel_k.cc)
target_include_directories(check_sph_bessel_k PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_bessel_k quadmath)
run_check(check_sph_bessel_k)

add_executable(check_sph_hankel_1 check/check_sph_hankel_1.cc)
target_include_directories(check_sph_hankel_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_hankel_1 quadmath)
run_check(check_sph_hankel_1)

add_executable(check_sph_hankel_2 check/check_sph_hankel_2.cc)
target_include_directories(check_sph_hankel_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_hankel_2 quadmath)
run_check(check_sph_hankel_2)

add_executable(check_sph_harmonic check/check_sph_harmonic.cc)
target_include_directories(check_sph_harmonic PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_harmonic quadmath)
run_check(check_sph_harmonic)

add_executable(check_sph_legendre check/check_sph_legendre.cc)
target_include_directories(check_sph_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_legendre quadmath)
run_check(check_sph_legendre)

add_executable(check_sph_neumann check/check_sph_neumann.cc)
target_include_directories(check_sph_neumann PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_sph_neumann quadmath)
run_check(check_sph_neumann)

add_executable(check_stirling_1 check/check_stirling_1.cc)
target_include_directories(check_stirling_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_stirling_1 quadmath)
run_check(check_stirling_1)

add_executable(check_stirling_2 check/check_stirling_2.cc)
target_include_directories(check_stirling_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_stirling_2 quadmath)
run_check(check_stirling_2)

add_executable(check_tgamma check/check_tgamma.cc)
target_include_directories(check_tgamma PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tgamma quadmath)
run_check(check_tgamma)

add_executable(check_tgamma_lower check/check_tgamma_lower.cc)
target_include_directories(check_tgamma_lower PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tgamma_lower quadmath)
run_check(check_tgamma_lower)

add_executable(check_theta_1 check/check_theta_1.cc)
target_include_directories(check_theta_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_theta_1 quadmath)
run_check(check_theta_1)

add_executable(check_theta_2 check/check_theta_2.cc)
target_include_directories(check_theta_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_theta_2 quadmath)
run_check(check_theta_2)

add_executable(check_theta_3 check/check_theta_3.cc)
target_include_directories(check_theta_3 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_theta_3 quadmath)
run_check(check_theta_3)

add_executable(check_theta_4 check/check_theta_4.cc)
target_include_directories(check_theta_4 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_theta_4 quadmath)
run_check(check_theta_4)

add_executable(check_tr1_assoc_laguerre check/check_tr1_assoc_laguerre.cc)
target_include_directories(check_tr1_assoc_laguerre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_assoc_laguerre quadmath)
run_check(check_tr1_assoc_laguerre)

add_executable(check_tr1_assoc_legendre check/check_tr1_assoc_legendre.cc)
target_include_directories(check_tr1_assoc_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_assoc_legendre quadmath)
run_check(check_tr1_assoc_legendre)

add_executable(check_tr1_beta check/check_tr1_beta.cc)
target_include_directories(check_tr1_beta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_beta quadmath)
run_check(check_tr1_beta)

add_executable(check_tr1_comp_ellint_1 check/check_tr1_comp_ellint_1.cc)
target_include_directories(check_tr1_comp_ellint_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_comp_ellint_1 quadmath)
run_check(check_tr1_comp_ellint_1)

add_executable(check_tr1_comp_ellint_2 check/check_tr1_comp_ellint_2.cc)
target_include_directories(check_tr1_comp_ellint_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_comp_ellint_2 quadmath)
run_check(check_tr1_comp_ellint_2)

add_executable(check_tr1_comp_ellint_3 check/check_tr1_comp_ellint_3.cc)
target_include_directories(check_tr1_comp_ellint_3 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_comp_ellint_3 quadmath)
run_check(check_tr1_comp_ellint_3)

add_executable(check_tr1_conf_hyperg check/check_tr1_conf_hyperg.cc)
target_include_directories(check_tr1_conf_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_conf_hyperg quadmath)
run_check(check_tr1_conf_hyperg)

add_executable(check_tr1_cyl_bessel_i check/check_tr1_cyl_bessel_i.cc)
target_include_directories(check_tr1_cyl_bessel_i PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_cyl_bessel_i quadmath)
run_check(check_tr1_cyl_bessel_i)

add_executable(check_tr1_cyl_bessel_j check/check_tr1_cyl_bessel_j.cc)
target_include_directories(check_tr1_cyl_bessel_j PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_cyl_bessel_j quadmath)
run_check(check_tr1_cyl_bessel_j)

add_executable(check_tr1_cyl_bessel_k check/check_tr1_cyl_bessel_k.cc)
target_include_directories(check_tr1_cyl_bessel_k PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_cyl_bessel_k quadmath)
run_check(check_tr1_cyl_bessel_k)

add_executable(check_tr1_cyl_neumann check/check_tr1_cyl_neumann.cc)
target_include_directories(check_tr1_cyl_neumann PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_cyl_neumann quadmath)
run_check(check_tr1_cyl_neumann)

add_executable(check_tr1_ellint_1 check/check_tr1_ellint_1.cc)
target_include_directories(check_tr1_ellint_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_ellint_1 quadmath)
run_check(check_tr1_ellint_1)

add_executable(check_tr1_ellint_2 check/check_tr1_ellint_2.cc)
target_include_directories(check_tr1_ellint_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_ellint_2 quadmath)
run_check(check_tr1_ellint_2)

add_executable(check_tr1_ellint_3 check/check_tr1_ellint_3.cc)
target_include_directories(check_tr1_ellint_3 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_ellint_3 quadmath)
run_check(check_tr1_ellint_3)

add_executable(check_tr1_expint check/check_tr1_expint.cc)
target_include_directories(check_tr1_expint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_expint quadmath)
run_check(check_tr1_expint)

add_executable(check_tr1_hermite check/check_tr1_hermite.cc)
target_include_directories(check_tr1_hermite PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_hermite quadmath)
run_check(check_tr1_hermite)

add_executable(check_tr1_hyperg check/check_tr1_hyperg.cc)
target_include_directories(check_tr1_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_hyperg quadmath)
run_check(check_tr1_hyperg)

add_executable(check_tr1_laguerre check/check_tr1_laguerre.cc)
target_include_directories(check_tr1_laguerre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_laguerre quadmath)
run_check(check_tr1_laguerre)

add_executable(check_tr1_legendre check/check_tr1_legendre.cc)
target_include_directories(check_tr1_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_legendre quadmath)
run_check(check_tr1_legendre)

add_executable(check_tr1_riemann_zeta check/check_tr1_riemann_zeta.cc)
target_include_directories(check_tr1_riemann_zeta PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_riemann_zeta quadmath)
run_check(check_tr1_riemann_zeta)

add_executable(check_tr1_sph_bessel check/check_tr1_sph_bessel.cc)
target_include_directories(check_tr1_sph_bessel PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_sph_bessel quadmath)
run_check(check_tr1_sph_bessel)

add_executable(check_tr1_sph_legendre check/check_tr1_sph_legendre.cc)
target_include_directories(check_tr1_sph_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_sph_legendre quadmath)
run_check(check_tr1_sph_legendre)

add_executable(check_tr1_sph_neumann check/check_tr1_sph_neumann.cc)
target_include_directories(check_tr1_sph_neumann PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_tr1_sph_neumann quadmath)
run_check(check_tr1_sph_neumann)

add_executable(check_zernike check/check_zernike.cc)
target_include_directories(check_zernike PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(check_zernike quadmath)
run_check(check_zernike)

add_executable(complex_airy_ai check/complex_airy_ai.cc)
target_include_directories(complex_airy_ai PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_airy_ai quadmath)
run_check(complex_airy_ai)

add_executable(complex_airy_bi check/complex_airy_bi.cc)
target_include_directories(complex_airy_bi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_airy_bi quadmath)
run_check(complex_airy_bi)

add_executable(complex_ellint_rc check/complex_ellint_rc.cc)
target_include_directories(complex_ellint_rc PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_ellint_rc quadmath)
run_check(complex_ellint_rc)

add_executable(complex_ellint_rd check/complex_ellint_rd.cc)
target_include_directories(complex_ellint_rd PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_ellint_rd quadmath)
run_check(complex_ellint_rd)

add_executable(complex_ellint_rf check/complex_ellint_rf.cc)
target_include_directories(complex_ellint_rf PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_ellint_rf quadmath)
run_check(complex_ellint_rf)

add_executable(complex_ellint_rg check/complex_ellint_rg.cc)
target_include_directories(complex_ellint_rg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_ellint_rg quadmath)
run_check(complex_ellint_rg)

add_executable(complex_ellint_rj check/complex_ellint_rj.cc)
target_include_directories(complex_ellint_rj PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(complex_ellint_rj quadmath)
run_check(complex_ellint_rj)

add_executable(deathmatch_comp_ellint check/deathmatch_comp_ellint.cc)
target_include_directories(deathmatch_comp_ellint PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(deathmatch_comp_ellint quadmath)

add_executable(deathmatch_conf_hyperg check/deathmatch_conf_hyperg.cc)
target_include_directories(deathmatch_conf_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(deathmatch_conf_hyperg quadmath)

add_executable(deathmatch_conf_hyperg_lim check/deathmatch_conf_hyperg_lim.cc)
target_include_directories(deathmatch_conf_hyperg_lim PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(deathmatch_conf_hyperg_lim quadmath)
run_check(deathmatch_conf_hyperg_lim)

add_executable(deathmatch_hyperg check/deathmatch_hyperg.cc)
target_include_directories(deathmatch_hyperg PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(deathmatch_hyperg quadmath)
run_check(deathmatch_hyperg)

add_executable(origin_cyl_bessel_j check/origin_cyl_bessel_j.cc)
target_include_directories(origin_cyl_bessel_j PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(origin_cyl_bessel_j quadmath)
run_check(origin_cyl_bessel_j)

add_executable(origin_cyl_neumann check/origin_cyl_neumann.cc)
target_include_directories(origin_cyl_neumann PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(origin_cyl_neumann quadmath)
run_check(origin_cyl_neumann)

add_executable(pr56216_cyl_bessel_i check/pr56216_cyl_bessel_i.cc)
target_include_directories(pr56216_cyl_bessel_i PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr56216_cyl_bessel_i quadmath)
run_check(pr56216_cyl_bessel_i)

add_executable(pr56216_cyl_hankel_1 check/pr56216_cyl_hankel_1.cc)
target_include_directories(pr56216_cyl_hankel_1 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr56216_cyl_hankel_1 quadmath)
run_check(pr56216_cyl_hankel_1)

add_executable(pr56216_cyl_hankel_2 check/pr56216_cyl_hankel_2.cc)
target_include_directories(pr56216_cyl_hankel_2 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr56216_cyl_hankel_2 quadmath)
run_check(pr56216_cyl_hankel_2)

add_executable(pr68397 check/pr68397.cc)
target_include_directories(pr68397 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr68397 quadmath)
run_check(pr68397)

add_executable(pr83140 check/pr83140.cc)
target_include_directories(pr83140 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr83140 quadmath)
run_check(pr83140)

add_executable(pr86655_assoc_legendre check/pr86655_assoc_legendre.cc)
target_include_directories(pr86655_assoc_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr86655_assoc_legendre quadmath)
run_check(pr86655_assoc_legendre)

add_executable(pr86655_sph_legendre check/pr86655_sph_legendre.cc)
target_include_directories(pr86655_sph_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr86655_sph_legendre quadmath)
run_check(pr86655_sph_legendre)

add_executable(pr86655_tr1_assoc_legendre check/pr86655_tr1_assoc_legendre.cc)
target_include_directories(pr86655_tr1_assoc_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr86655_tr1_assoc_legendre quadmath)
run_check(pr86655_tr1_assoc_legendre)

add_executable(pr86655_tr1_sph_legendre check/pr86655_tr1_sph_legendre.cc)
target_include_directories(pr86655_tr1_sph_legendre PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(pr86655_tr1_sph_legendre quadmath)
run_check(pr86655_tr1_sph_legendre)

# This needs C++14!
#add_executable(tr1_airy check/tr1_airy.cc)
#target_include_directories(tr1_airy PRIVATE
#  .
#  ${TR29124_INCLUDE_DIRS}
#  ../wrappers)
#target_link_libraries(tr1_airy quadmath)
#run_check(tr1_airy)

add_executable(tr1_pr68397 check/tr1_pr68397.cc)
target_include_directories(tr1_pr68397 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(tr1_pr68397 quadmath)
run_check(tr1_pr68397)

add_executable(tr1_pr83140 check/tr1_pr83140.cc)
target_include_directories(tr1_pr83140 PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../wrappers)
target_link_libraries(tr1_pr83140 quadmath)
run_check(tr1_pr83140)

