cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_bernoulli_functions
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

add_executable(build_bernoulli_2n_table build_bernoulli_2n_table.cpp)
target_include_directories(build_bernoulli_2n_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(build_bernoulli_2n_table quadmath mpfr)
add_test(NAME run_build_bernoulli_2n_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_bernoulli_2n_table > output/build_bernoulli_2n_table.txt")

add_executable(test_bernoulli test_bernoulli.cpp)
target_include_directories(test_bernoulli PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_float128/include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_bernoulli quadmath
  -L../../wrappers -lwrap_burkhardt
  -L/home/ed/bin/lib64 -lgfortran)
add_test(NAME run_test_bernoulli COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bernoulli > output/test_bernoulli.txt")
