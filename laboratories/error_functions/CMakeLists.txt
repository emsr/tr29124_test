cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_error_functions
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

add_executable(build_inv_erf_coefs build_inv_erf_coefs.cpp)
target_include_directories(build_inv_erf_coefs PRIVATE
  $ENV{HOME}/mpreal
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_inv_erf_coefs
  quadmath
  gmp mpfr)
add_test(NAME run_build_inv_erf_coefs COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_inv_erf_coefs > build_inv_erf_coefs.txt")

add_executable(test_dawson test_dawson.cpp)
target_include_directories(test_dawson PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_dawson
  quadmath
  -L../../wrappers -lwrap_boost)
add_test(NAME run_test_dawson COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dawson > test_dawson.txt")

add_executable(test_erfc test_erfc.cpp)
target_include_directories(test_erfc PRIVATE
  ../../cxx_continued_fractions/include
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_erfc quadmath)
add_test(NAME run_test_erfc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_erfc > test_erfc.txt")

add_executable(test_experfc test_experfc.cpp)
target_include_directories(test_experfc PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_continued_fractions/include
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_experfc
  quadmath
  mpfr gmp)
add_test(NAME run_test_experfc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_experfc > test_experfc.txt")

add_executable(test_faddeeva test_faddeeva.cpp)
target_include_directories(test_faddeeva PRIVATE
  .
  ../../wrappers)
add_test(NAME run_test_faddeeva COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_faddeeva > test_faddeeva.txt")

add_executable(test_fresnel test_fresnel.cpp)
target_include_directories(test_fresnel PRIVATE
  .
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_fresnel quadmath)
add_test(NAME run_test_fresnel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_fresnel > test_fresnel.txt")

add_executable(test_inv_erf test_inv_erf.cpp)
target_compile_definitions(test_inv_erf PRIVATE STANDALONE)
target_include_directories(test_inv_erf PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
)
target_link_libraries(test_inv_erf quadmath)
add_test(NAME run_test_inv_erf COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_erf > test_inv_erf.txt")

add_executable(test_owens_t test_owens_t.cpp)
target_include_directories(test_owens_t PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
)
target_link_libraries(test_owens_t
  quadmath
  -L../../wrappers -lwrap_boost
)
add_test(NAME run_test_owens_t COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_owens_t > test_owens_t.txt")
