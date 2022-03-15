cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_elementary_functions
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

add_custom_target(make_lab_elementary_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

if (HAVE_MPREAL)
  add_executable(build_atan_table build_atan_table.cpp)
  target_link_libraries(build_atan_table ${EMSR_LIBRARIES} mpreal gmp mpfr quadmath)
  add_test(NAME run_build_atan_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_atan_table > output/build_atan_table.txt")

  add_executable(build_cordic build_cordic.cpp)
  target_link_libraries(build_cordic ${EMSR_LIBRARIES} cxx_float128 mpreal gmp mpfr quadmath)
  add_test(NAME run_build_cordic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_cordic > output/build_cordic.txt")

  add_executable(build_log_table build_log_table.cpp)
  target_link_libraries(build_log_table ${EMSR_LIBRARIES} mpreal gmp mpfr quadmath)
  add_test(NAME run_build_log_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_log_table > output/build_log_table.txt")

  add_executable(build_sincos_tables build_sincos_tables.cpp)
  target_link_libraries(build_sincos_tables ${EMSR_LIBRARIES} mpreal gmp mpfr quadmath)
  add_test(NAME run_build_sincos_tables COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_sincos_tables > output/build_sincos_tables.txt")

  add_executable(build_sqrt_table build_sqrt_table.cpp)
  target_link_libraries(build_sqrt_table ${EMSR_LIBRARIES} mpreal gmp mpfr quadmath)
  add_test(NAME run_build_sqrt_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_sqrt_table > output/build_sqrt_table.txt")
endif(HAVE_MPREAL)

add_executable(test_cordic test_cordic.cpp)
target_include_directories(test_cordic PRIVATE
  .)
add_test(NAME run_test_cordic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cordic > output/test_cordic.txt")

add_executable(test_exp test_exp.cpp)
target_link_libraries(test_exp ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_exp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_exp > output/test_exp.txt")

add_executable(test_gudermannian test_gudermannian.cpp)
target_link_libraries(test_gudermannian ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_gudermannian COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gudermannian > output/test_gudermannian.txt")

add_executable(test_lambert_w test_lambert_w.cpp)
target_link_libraries(test_lambert_w ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_lambert_w COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lambert_w > output/test_lambert_w.txt")

add_executable(test_log test_log.cpp)
add_test(NAME run_test_log COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_log > output/test_log.txt")

add_executable(test_notsospecfun test_notsospecfun.cpp)
target_link_libraries(test_notsospecfun ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_notsospecfun COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_notsospecfun > output/test_notsospecfun.txt")

add_executable(test_pow test_pow.cpp)
add_test(NAME run_test_pow COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_pow > output/test_pow.txt")

add_executable(test_reperiodized_hyper test_reperiodized_hyper.cpp)
target_link_libraries(test_reperiodized_hyper ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_reperiodized_hyper COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_reperiodized_hyper > output/test_reperiodized_hyper.txt")

if (HAVE_BOOST)
  add_executable(test_reperiodized_trig test_reperiodized_trig.cpp)
  target_link_libraries(test_reperiodized_trig ${EMSR_LIBRARIES} cxx_float128 quadmath wrap_boost)
  add_test(NAME run_test_reperiodized_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_reperiodized_trig > output/test_reperiodized_trig.txt")

  add_executable(test_sinus_cardinal test_sinus_cardinal.cpp)
  target_link_libraries(test_sinus_cardinal ${EMSR_LIBRARIES} cxx_float128 quadmath wrap_gsl wrap_boost)
  add_test(NAME run_test_sinus_cardinal COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sinus_cardinal > output/test_sinus_cardinal.txt")
endif(HAVE_BOOST)

add_executable(test_sincos test_sincos.cpp)
target_link_libraries(test_sincos ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_sincos COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sincos > output/test_sincos.txt")

add_executable(test_sqrt test_sqrt.cpp)
target_link_libraries(test_sqrt ${EMSR_LIBRARIES})
add_test(NAME run_test_sqrt COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sqrt > output/test_sqrt.txt")

add_executable(test_static_root test_static_root.cpp)
add_test(NAME run_test_static_root COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_static_root > output/test_static_root.txt")

add_executable(test_trig test_trig.cpp)
target_link_libraries(test_trig ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_trig > output/test_trig.txt")

add_executable(test_wright_omega test_wright_omega.cpp)
target_link_libraries(test_wright_omega ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_wright_omega COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_wright_omega > output/test_wright_omega.txt")

add_executable(test_zeta_trig test_zeta_trig.cpp)
target_link_libraries(test_zeta_trig ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_zeta_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_zeta_trig > output/test_zeta_trig.txt")
