cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_elementary_functions
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

add_executable(build_atan_table build_atan_table.cpp)
target_include_directories(build_atan_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_atan_table gmp mpfr quadmath)
add_test(NAME run_build_atan_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_atan_table > output/build_atan_table.txt")

add_executable(build_cordic build_cordic.cpp)
target_include_directories(build_cordic PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_cordic gmp mpfr quadmath)
add_test(NAME run_build_cordic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_cordic > output/build_cordic.txt")

add_executable(build_log_table build_log_table.cpp)
target_include_directories(build_log_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_log_table gmp mpfr quadmath)
add_test(NAME run_build_log_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_log_table > output/build_log_table.txt")

add_executable(build_sincos_tables build_sincos_tables.cpp)
target_include_directories(build_sincos_tables PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_sincos_tables gmp mpfr quadmath)
add_test(NAME run_build_sincos_tables COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_sincos_tables > output/build_sincos_tables.txt")

add_executable(build_sqrt_table build_sqrt_table.cpp)
target_include_directories(build_sqrt_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_sqrt_table gmp mpfr quadmath)
add_test(NAME run_build_sqrt_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_sqrt_table > output/build_sqrt_table.txt")

add_executable(test_cordic test_cordic.cpp)
target_include_directories(test_cordic PRIVATE
  .)
add_test(NAME run_test_cordic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cordic > output/test_cordic.txt")

add_executable(test_exp test_exp.cpp)
target_include_directories(test_exp PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_exp quadmath)
add_test(NAME run_test_exp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_exp > output/test_exp.txt")

add_executable(test_gudermannian test_gudermannian.cpp)
target_include_directories(test_gudermannian PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_gudermannian quadmath)
add_test(NAME run_test_gudermannian COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gudermannian > output/test_gudermannian.txt")

add_executable(test_lambert_w test_lambert_w.cpp)
target_include_directories(test_lambert_w PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_lambert_w quadmath)
add_test(NAME run_test_lambert_w COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lambert_w > output/test_lambert_w.txt")

add_executable(test_log test_log.cpp)
add_test(NAME run_test_log COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_log > output/test_log.txt")

add_executable(test_notsospecfun test_notsospecfun.cpp)
target_include_directories(test_notsospecfun PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_notsospecfun quadmath)
add_test(NAME run_test_notsospecfun COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_notsospecfun > output/test_notsospecfun.txt")

add_executable(test_pow test_pow.cpp)
add_test(NAME run_test_pow COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_pow > output/test_pow.txt")

add_executable(test_reperiodized_hyper test_reperiodized_hyper.cpp)
target_include_directories(test_reperiodized_hyper PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_reperiodized_hyper quadmath)
add_test(NAME run_test_reperiodized_hyper COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_reperiodized_hyper > output/test_reperiodized_hyper.txt")

add_executable(test_reperiodized_trig test_reperiodized_trig.cpp)
target_include_directories(test_reperiodized_trig PRIVATE
  ../../wrappers
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_reperiodized_trig quadmath -L../../wrappers -lwrap_boost)
add_test(NAME run_test_reperiodized_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_reperiodized_trig > output/test_reperiodized_trig.txt")

add_executable(test_sincos test_sincos.cpp)
target_include_directories(test_sincos PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_sincos quadmath)
add_test(NAME run_test_sincos COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sincos > output/test_sincos.txt")

add_executable(test_sinus_cardinal test_sinus_cardinal.cpp)
target_include_directories(test_sinus_cardinal PRIVATE
  ../../wrappers
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_sinus_cardinal
   quadmath
   -L../../wrappers -lwrap_gsl -lwrap_boost)
add_test(NAME run_test_sinus_cardinal COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sinus_cardinal > output/test_sinus_cardinal.txt")

add_executable(test_sqrt test_sqrt.cpp)
target_include_directories(test_sqrt PRIVATE
  ../../cxx_math_constants/include)
add_test(NAME run_test_sqrt COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sqrt > output/test_sqrt.txt")

add_executable(test_static_root test_static_root.cpp)
add_test(NAME run_test_static_root COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_static_root > output/test_static_root.txt")

add_executable(test_trig test_trig.cpp)
target_include_directories(test_trig PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_trig quadmath)
add_test(NAME run_test_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_trig > output/test_trig.txt")

add_executable(test_wright_omega test_wright_omega.cpp)
target_include_directories(test_wright_omega PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_wright_omega quadmath)
add_test(NAME run_test_wright_omega COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_wright_omega > output/test_wright_omega.txt")

add_executable(test_zeta_trig test_zeta_trig.cpp)
target_include_directories(test_zeta_trig PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_zeta_trig quadmath)
add_test(NAME run_test_zeta_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_zeta_trig > output/test_zeta_trig.txt")
