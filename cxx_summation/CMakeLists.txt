cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_summation
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_summation_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

# Header-only library.
add_library(cxx_summation INTERFACE)
target_include_directories(cxx_summation INTERFACE include)

add_executable(test_summation test_summation.cpp)
target_include_directories(test_summation PRIVATE
  ../include)
target_link_libraries(test_summation
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_polynomial
  cxx_continued_fractions
  cxx_summation
  cxx_integration
  cxx_float128
  quadmath)
