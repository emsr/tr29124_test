cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_integer_functions
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

add_executable(test_factorization test_factorization.cpp)
add_test(NAME run_test_factorization COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_factorization > output/test_factorization.txt")

add_executable(test_prime_factors test_prime_factors.cpp)
target_include_directories(test_prime_factors PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_continued_fractions/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_prime_factors quadmath)
add_test(NAME run_test_prime_factors COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_prime_factors > output/test_prime_factors.txt")
