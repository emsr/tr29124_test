cmake_minimum_required (VERSION 3.10)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX
)

enable_testing()

set(EMSR_LIBRARIES
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_continued_fractions
  cxx_summation
  cxx_polynomial
  cxx_integration
  cxx_special_functions)

add_subdirectory(airy_functions)
add_subdirectory(appell_functions)
add_subdirectory(basic_functions)
add_subdirectory(bernoulli_functions)
add_subdirectory(bessel_functions)
add_subdirectory(beta_functions)
add_subdirectory(complex_tools)
add_subdirectory(constants)
add_subdirectory(coulomb_functions)
add_subdirectory(distribution_functions)
add_subdirectory(elementary_functions)
add_subdirectory(elliptic_integrals)
add_subdirectory(error_functions)
add_subdirectory(exponential_integrals)
add_subdirectory(floating_point_tools)
add_subdirectory(gamma_functions)
#add_subdirectory(heun_functions)
add_subdirectory(hypergeometric_functions)
add_subdirectory(integer_functions)
add_subdirectory(mathieu_functions)
add_subdirectory(mittag_leffler_functions)
add_subdirectory(norm_functions)
add_subdirectory(orthogonal_polynomials)
add_subdirectory(parabolic_cylinder_functions)
add_subdirectory(theta_functions)
add_subdirectory(zeta_functions)
add_subdirectory(recursive_functions)
add_subdirectory(multiprecision)
add_subdirectory(harmonics)
add_subdirectory(polynomials)
