cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_beta_functions
  VERSION 0.5
  LANGUAGES CXX Fortran
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

add_executable(test_beta test_beta.cpp)
target_include_directories(test_beta PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_beta quadmath -L../../wrappers -lwrap_boost)
add_test(NAME run_test_beta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_beta > output/test_beta.txt")

add_executable(test_beta_inc test_beta_inc.cpp)
target_include_directories(test_beta_inc PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_beta_inc quadmath)
add_test(NAME run_test_beta_inc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_beta_inc > output/test_beta_inc.txt")

add_executable(test_inv_ibeta test_inv_ibeta.cpp)
target_include_directories(test_inv_ibeta PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_root_search/include)
target_link_libraries(test_inv_ibeta quadmath)
add_test(NAME run_test_inv_ibeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_ibeta > output/test_inv_ibeta.txt")
