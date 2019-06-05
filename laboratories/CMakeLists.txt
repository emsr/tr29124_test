cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_CXX_FLAGS -Wno-psabi)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_subdirectory(airy_functions)
add_subdirectory(appell_functions)
add_subdirectory(basic_functions)
add_subdirectory(bernoulli_functions)
add_subdirectory(bessel_functions)
add_subdirectory(beta_functions)
add_subdirectory(complex_tools)
add_subdirectory(constants)
add_subdirectory(coulomb_functions)
add_subdirectory(distributions)
add_subdirectory(elementary_functions)
add_subdirectory(elliptic_integrals)
add_subdirectory(error_functions)
add_subdirectory(exponential_integrals)
add_subdirectory(floating_point_tools)
add_subdirectory(gamma_functions)
add_subdirectory(heun_functions)
add_subdirectory(hypergeometric_functions)
add_subdirectory(integer_functions)
add_subdirectory(mathieu_functions)
add_subdirectory(mittag_leffler_functions)
add_subdirectory(norm_functions)
add_subdirectory(orthogonal_polynomials)
add_subdirectory(parabolic_cylinder_functions)
add_subdirectory(theta_functions)
add_subdirectory(zeta_functions)
