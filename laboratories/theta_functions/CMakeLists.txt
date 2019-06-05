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

add_executable(test_jacobi_ellint test_jacobi_ellint.cpp)
target_include_directories(test_jacobi_ellint PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_jacobi_ellint quadmath -L../../wrappers/debug -lwrap_gsl)

add_executable(test_jacobi_inv test_jacobi_inv.cpp)
target_include_directories(test_jacobi_inv PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_jacobi_inv quadmath)

add_executable(test_jacobi_theta test_jacobi_theta.cpp)
target_include_directories(test_jacobi_theta PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_jacobi_theta quadmath -L../../wrappers/debug -lwrap_gsl)

add_executable(test_theta test_theta.cpp)
target_include_directories(test_theta PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_theta quadmath)

add_executable(test_weierstrass_ellint test_weierstrass_ellint.cpp)
target_include_directories(test_weierstrass_ellint PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_weierstrass_ellint quadmath)
