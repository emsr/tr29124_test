cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX Fortran
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_FORTRAN_COMPILER "/home/ed/bin/bin/gfortran")
set(CMAKE_CXX_FLAGS -Wno-psabi)

enable_testing()

add_executable(run_coulfg run_coulfg.cpp coulfg.cpp)
#target_include_directories(run_coulfg PRIVATE )
#target_link_libraries(run_coulfg )

add_executable(test_coulomb test_coulomb.cpp)
#target_include_directories(test_coulomb PRIVATE )
#target_link_libraries(test_coulomb )

add_executable(test_hydrogen test_hydrogen.cpp)
target_include_directories(test_hydrogen PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_hydrogen quadmath)
