cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_norm_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_executable(test_abs test_abs.cpp)
add_test(NAME run_test_abs COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_abs > output/test_abs.txt")

add_executable(test_arith_geom_mean test_arith_geom_mean.cpp)
target_include_directories(test_arith_geom_mean PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_float128/include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../cxx_complex_math/include
  ../../wrappers)
target_link_libraries(test_arith_geom_mean
  quadmath
  -L../../wrappers -lwrap_boost)
add_test(NAME run_test_arith_geom_mean COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_arith_geom_mean > output/test_arith_geom_mean.txt")

add_executable(test_hypot test_hypot.cpp)
target_include_directories(test_hypot PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_continued_fractions/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_hypot quadmath)
add_test(NAME run_test_hypot COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hypot > output/test_hypot.txt")

add_executable(test_logsumexp test_logsumexp.cpp)
add_test(NAME run_test_logsumexp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_logsumexp > output/test_logsumexp.txt")

add_executable(test_power_mean test_power_mean.cpp)
add_test(NAME run_test_power_mean COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_power_mean > output/test_power_mean.txt")

add_executable(test_power_norm test_power_norm.cpp)
add_test(NAME run_test_power_norm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_power_norm > output/test_power_norm.txt")
