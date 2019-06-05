cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_norm_functions
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

add_executable(test_abs test_abs.cpp)
#target_include_directories(test_abs PRIVATE )
#target_link_libraries(test_abs )

add_executable(test_arith_geom_mean test_arith_geom_mean.cpp)
target_include_directories(test_arith_geom_mean PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../wrappers)
target_link_libraries(test_arith_geom_mean
  quadmath
  -L../../wrappers/debug -lwrap_boost)

add_executable(test_hypot test_hypot.cpp)
target_include_directories(test_hypot PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_hypot quadmath)

add_executable(test_logsumexp test_logsumexp.cpp)
#target_include_directories(test_logsumexp PRIVATE )
#target_link_libraries(test_logsumexp )

add_executable(test_power_mean test_power_mean.cpp)
#target_include_directories(test_power_mean PRIVATE )
#target_link_libraries(test_power_mean )

add_executable(test_power_norm test_power_norm.cpp)
#target_include_directories(test_power_norm PRIVATE )
#target_link_libraries(test_power_norm )
