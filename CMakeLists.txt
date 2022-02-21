cmake_minimum_required (VERSION 3.10)

project(
  cxx_maths
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

enable_testing()

# See if we can have Boost...

#find_package(Boost REQUIRED)
find_package(Boost)

set(HAVE_BOOST FALSE CACHE BOOL "Have working boost" FORCE)
if (Boost_FOUND)
  set(HAVE_BOOST TRUE CACHE BOOL "Have working boost" FORCE)
  message(STATUS "Boost found!")
else()
  message(WARNING "Boost not found!")
endif(Boost_FOUND)

# See if we can have mpreal...

find_package(GMP)
if (GMP_FOUND)
  message(STATUS "GMP found!")
else()
  message(WARNING "GMP not found!")
endif()

find_package(MPFR)
if (MPFR_FOUND)
  message(STATUS "MPFR found!")
else()
  message(WARNING "MPFR not found!")
endif()

set(HAVE_MPREAL FALSE CACHE BOOL "Have working mpreal" FORCE)
if (GMP_FOUND AND MPFR_FOUND)
  set(HAVE_MPREAL TRUE CACHE BOOL "Have working mpreal" FORCE)
  message(STATUS "MPReal enabled!")
else()
  message(WARNING "MPReal not enabled!")
endif()

# Add all the subdirectories...

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
add_subdirectory(cxx_matrix_math)

add_subdirectory(3rdparty)

# Actual submodules
add_subdirectory(cxx_math_constants)
add_subdirectory(cxx_polynomial)
add_subdirectory(cxx_integration)
