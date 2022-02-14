cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_gamma_functions
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

add_custom_target(make_lab_gamma_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(build_gamma_lanczos build_gamma_lanczos.cpp)
target_include_directories(build_gamma_lanczos PRIVATE
  ../../cxx_float128/include
)
target_link_libraries(build_gamma_lanczos ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_build_gamma_lanczos COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_gamma_lanczos > output/build_gamma_lanczos.txt")

if (HAVE_MPREAL)
  add_executable(build_gamma_recip build_gamma_recip.cpp)
  target_link_libraries(build_gamma_recip ${EMSR_LIBRARIES} quadmath mpreal mpfr gmp)
  add_test(NAME run_build_gamma_recip COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_gamma_recip > output/build_gamma_recip.txt")
endif(HAVE_MPREAL)

add_executable(build_gamma_spouge build_gamma_spouge.cpp)
target_link_libraries(build_gamma_spouge ${EMSR_LIBRARIES} cxx_float128 quadmath)
#FIXME add_test(NAME run_build_gamma_spouge COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_gamma_spouge > output/build_gamma_spouge.txt")

add_executable(plot_gamma plot_gamma.cpp)
target_link_libraries(plot_gamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_plot_gamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/plot_gamma > output/plot_gamma.txt")

add_executable(test_binet test_binet.cpp)
target_link_libraries(test_binet ${EMSR_LIBRARIES} cxx_rational quadmath)
add_test(NAME run_test_binet COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_binet > output/test_binet.txt")

add_executable(test_binet_float test_binet_float.cpp)
target_link_libraries(test_binet_float ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_binet_float COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_binet_float > output/test_binet_float.txt")

add_executable(test_complex_gamma test_complex_gamma.cpp)
target_link_libraries(test_complex_gamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_complex_gamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_complex_gamma > output/test_complex_gamma.txt")

add_executable(test_digamma test_digamma.cpp)
target_link_libraries(test_digamma
  ${EMSR_LIBRARIES}
  quadmath wrap_gsl)
add_test(NAME run_test_digamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_digamma > output/test_digamma.txt")

add_executable(test_factorial test_factorial.cpp)
target_link_libraries(test_factorial ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_factorial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_factorial > output/test_factorial.txt")

if (HAVE_BOOST)
  add_executable(test_falling_factorial test_falling_factorial.cpp)
  target_link_libraries(test_falling_factorial
    ${EMSR_LIBRARIES}
    quadmath wrap_boost)
  add_test(NAME run_test_falling_factorial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_falling_factorial > output/test_falling_factorial.txt")

  add_executable(test_gamma test_gamma.cpp)
  target_link_libraries(test_gamma quadmath ${EMSR_LIBRARIES} cxx_float128 wrap_boost)
  add_test(NAME run_test_gamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gamma > output/test_gamma.txt")

  add_executable(test_gamma_ratio test_gamma_ratio.cpp)
  target_link_libraries(test_gamma_ratio ${EMSR_LIBRARIES} cxx_float128 quadmath wrap_boost)
  add_test(NAME run_test_gamma_ratio COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gamma_ratio > output/test_gamma_ratio.txt")
endif(HAVE_BOOST)

add_executable(test_gamma_reciprocal test_gamma_reciprocal.cpp)
target_link_libraries(test_gamma_reciprocal ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_gamma_reciprocal COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gamma_reciprocal > output/test_gamma_reciprocal.txt")

add_executable(test_inv_gamma test_inv_gamma.cpp)
target_link_libraries(test_inv_gamma ${EMSR_LIBRARIES} cxx_float128 quadmath)
#FIXME add_test(NAME run_test_inv_gamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_gamma > output/test_inv_gamma.txt")

add_executable(test_inv_lgamma test_inv_lgamma.cpp)
target_include_directories(test_inv_lgamma PRIVATE
  .)
target_link_libraries(test_inv_lgamma ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_inv_lgamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_lgamma > output/test_inv_lgamma.txt")

if (HAVE_BOOST)
  add_executable(test_polygamma test_polygamma.cpp)
  target_link_libraries(test_polygamma ${EMSR_LIBRARIES} quadmath wrap_boost)
  add_test(NAME run_test_polygamma COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_polygamma > output/test_polygamma.txt")

  add_executable(test_rising_factorial test_rising_factorial.cpp)
  target_link_libraries(test_rising_factorial ${EMSR_LIBRARIES} quadmath wrap_boost)
  add_test(NAME run_test_rising_factorial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_rising_factorial > output/test_rising_factorial.txt")
endif(HAVE_BOOST)
