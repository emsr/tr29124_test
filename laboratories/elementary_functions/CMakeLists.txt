cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_elementary_functions
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

add_executable(build_atan_table build_atan_table.cpp)
target_include_directories(build_atan_table PRIVATE
  $ENV{HOME}/mpreal
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_atan_table gmp mpfr quadmath)

add_executable(build_cordic build_cordic.cpp)
target_include_directories(build_cordic PRIVATE
  $ENV{HOME}/mpreal
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_cordic gmp mpfr quadmath)

add_executable(build_log_table build_log_table.cpp)
target_include_directories(build_log_table PRIVATE
  $ENV{HOME}/mpreal
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_log_table gmp mpfr quadmath)

add_executable(build_sincos_tables build_sincos_tables.cpp)
target_include_directories(build_sincos_tables PRIVATE
  $ENV{HOME}/mpreal
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_sincos_tables gmp mpfr quadmath)

add_executable(build_sqrt_table build_sqrt_table.cpp)
target_include_directories(build_sqrt_table PRIVATE
  $ENV{HOME}/mpreal
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_sqrt_table gmp mpfr quadmath)

add_executable(test_cordic test_cordic.cpp)
target_include_directories(test_cordic PRIVATE
  .
)
#target_link_libraries(test_cordic )

add_executable(test_exp test_exp.cpp)
target_include_directories(test_exp PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_exp quadmath)

add_executable(test_gudermannian test_gudermannian.cpp)
target_include_directories(test_gudermannian PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_gudermannian quadmath)

add_executable(test_lambert_w test_lambert_w.cpp)
target_include_directories(test_lambert_w PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_lambert_w quadmath)

add_executable(test_log test_log.cpp)
#target_include_directories(test_log PRIVATE )
#target_link_libraries(test_log )

add_executable(test_notsospecfun test_notsospecfun.cpp)
target_include_directories(test_notsospecfun PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_notsospecfun quadmath)

add_executable(test_pow test_pow.cpp)
#target_include_directories(test_pow PRIVATE )
#target_link_libraries(test_pow )

add_executable(test_reperiodized_hyper test_reperiodized_hyper.cpp)
target_include_directories(test_reperiodized_hyper PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_reperiodized_hyper quadmath)

add_executable(test_reperiodized_trig test_reperiodized_trig.cpp)
target_include_directories(test_reperiodized_trig PRIVATE
  ../../wrappers
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_reperiodized_trig quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_sincos test_sincos.cpp)
target_include_directories(test_sincos PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_sincos quadmath)

add_executable(test_sinus_cardinal test_sinus_cardinal.cpp)
target_include_directories(test_sinus_cardinal PRIVATE
  ../../wrappers
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_sinus_cardinal
   quadmath
   -L../../wrappers/debug -lwrap_gsl -lwrap_boost)

add_executable(test_sqrt test_sqrt.cpp)
#target_include_directories(test_sqrt PRIVATE )
#target_link_libraries(test_sqrt )

add_executable(test_static_root test_static_root.cpp)
#target_include_directories(test_static_root PRIVATE )
#target_link_libraries(test_static_root )

add_executable(test_trig test_trig.cpp)
target_include_directories(test_trig PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_trig quadmath)

add_executable(test_wright_omega test_wright_omega.cpp)
target_include_directories(test_wright_omega PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_wright_omega quadmath)

add_executable(test_zeta_trig test_zeta_trig.cpp)
target_include_directories(test_zeta_trig PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_zeta_trig quadmath)
