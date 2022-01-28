cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_error_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_lab_error_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(build_inv_erf_coefs build_inv_erf_coefs.cpp)
target_link_libraries(build_inv_erf_coefs
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  mpreal gmp mpfr)
add_test(NAME run_build_inv_erf_coefs COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_inv_erf_coefs > output/build_inv_erf_coefs.txt")

add_executable(test_dawson test_dawson.cpp)
target_link_libraries(test_dawson
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  wrap_boost)
add_test(NAME run_test_dawson COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dawson > output/test_dawson.txt")

add_executable(test_erfc test_erfc.cpp)
target_link_libraries(test_erfc ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_erfc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_erfc > output/test_erfc.txt")

add_executable(test_experfc test_experfc.cpp)
target_link_libraries(test_experfc
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_test_experfc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_experfc > output/test_experfc.txt")

add_executable(test_faddeeva test_faddeeva.cpp)
target_link_libraries(test_faddeeva
  wrap_faddeeva)
add_test(NAME run_test_faddeeva COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_faddeeva > output/test_faddeeva.txt")

add_executable(test_fresnel test_fresnel.cpp)
target_include_directories(test_fresnel PRIVATE
  .)
target_link_libraries(test_fresnel
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath)
add_test(NAME run_test_fresnel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_fresnel > output/test_fresnel.txt")

add_executable(test_inv_erf test_inv_erf.cpp)
target_compile_definitions(test_inv_erf PRIVATE STANDALONE)
target_link_libraries(test_inv_erf
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath)
add_test(NAME run_test_inv_erf COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_erf > output/test_inv_erf.txt")

add_executable(test_owens_t test_owens_t.cpp)
target_link_libraries(test_owens_t
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  wrap_boost)
add_test(NAME run_test_owens_t COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_owens_t > output/test_owens_t.txt")
