cmake_minimum_required (VERSION 3.10)

project(
  cxx_maths
  VERSION 0.5
  LANGUAGES CXX
)

enable_testing()

add_subdirectory(wrappers)
add_subdirectory(libstdc++_support)
add_subdirectory(laboratories)

#add_subdirectory(cxx_chebyshev)# docs!
add_subdirectory(cxx_complex_math)
add_subdirectory(cxx_complex_utils)
add_subdirectory(cxx_continued_fractions)
#add_subdirectory(cxx_float128)
add_subdirectory(cxx_fp_utils)
add_subdirectory(cxx_interval)
#add_subdirectory(cxx_rational)# docs!
add_subdirectory(cxx_root_search)
add_subdirectory(cxx_summation)
add_subdirectory(cxx_traits_utils)

# Actual submodules
#add_subdirectory(cxx_math_constants)# docs!
add_subdirectory(polynomial)
add_subdirectory(quadrature)
