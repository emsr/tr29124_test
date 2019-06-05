cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_beta_functions
  VERSION 0.5
  LANGUAGES CXX Fortran
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_FORTRAN_COMPILER "/home/ed/bin/bin/gfortran")
set(CMAKE_CXX_FLAGS -Wno-psabi)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_executable(test_beta test_beta.cpp)
target_include_directories(test_beta PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_beta quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_beta_inc test_beta_inc.cpp)
target_include_directories(test_beta_inc PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_beta_inc quadmath)

add_executable(test_inv_ibeta test_inv_ibeta.cpp)
target_include_directories(test_inv_ibeta PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_root_search/include)
target_link_libraries(test_inv_ibeta quadmath)
