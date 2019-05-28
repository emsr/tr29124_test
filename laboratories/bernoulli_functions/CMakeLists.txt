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

enable_testing()

add_executable(build_bernoulli_2n_table build_bernoulli_2n_table.cpp)
target_include_directories(build_bernoulli_2n_table PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(build_bernoulli_2n_table quadmath mpfr)

add_executable(test_bernoulli test_bernoulli.cpp)
target_include_directories(test_bernoulli PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_bernoulli quadmath
  -L../../wrappers/debug -lwrap_burkhardt
  -L/home/ed/bin/lib64 -lgfortran)
