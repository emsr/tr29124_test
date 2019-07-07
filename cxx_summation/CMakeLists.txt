cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_summation
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_executable(test_summation test_summation.cpp)
target_include_directories(test_summation PRIVATE
  include
  ../include
  ../cxx_float128/include
  ../cxx_math_constants/include
  ../cxx_fp_utils/include
  ../cxx_complex_utils/include
  ../cxx_complex_math/include
  ../polynomial/include
  ../quadrature/include)
target_link_libraries(test_summation quadmath)
