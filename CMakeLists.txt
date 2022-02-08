cmake_minimum_required (VERSION 3.10)

project(
  cxx_maths
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

# Include BoostLib module .. Can't get it to work.
#set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/node_modules/boost-lib/cmake")
#include(BoostLib)
#require_boost_libs(">= 1.74.0" "")

enable_testing()

add_subdirectory(wrappers)
add_subdirectory(laboratories)

add_subdirectory(cxx_chebyshev)
add_subdirectory(cxx_complex_math)
add_subdirectory(cxx_complex_utils)
add_subdirectory(cxx_continued_fractions)
add_subdirectory(cxx_float128)
add_subdirectory(cxx_fp_utils)
add_subdirectory(cxx_interval)
add_subdirectory(cxx_numeric_limits)
add_subdirectory(cxx_rational)
add_subdirectory(cxx_root_search)
add_subdirectory(cxx_summation)
add_subdirectory(cxx_traits_utils)
add_subdirectory(cxx_special_functions)

add_subdirectory(3rdparty)

# Actual submodules
add_subdirectory(cxx_math_constants)
add_subdirectory(cxx_polynomial)
add_subdirectory(cxx_integration)
