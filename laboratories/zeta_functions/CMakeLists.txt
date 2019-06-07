cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_zeta_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_executable(build_etam1_table build_etam1_table.cpp)
target_include_directories(build_etam1_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_etam1_table
  quadmath
  mpfr gmp)

add_executable(build_nfact_zetanp1 build_nfact_zetanp1.cpp)
target_include_directories(build_nfact_zetanp1 PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_nfact_zetanp1
  quadmath
  mpfr gmp)

add_executable(build_zeta_deriv_table build_zeta_deriv_table.cpp)
target_include_directories(build_zeta_deriv_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_zeta_deriv_table
  quadmath
  mpfr gmp)

add_executable(build_zetahalfm1_table build_zetahalfm1_table.cpp)
target_include_directories(build_zetahalfm1_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_zetahalfm1_table
  quadmath
  mpfr gmp)

add_executable(build_zetam1_table build_zetam1_table.cpp)
target_include_directories(build_zetam1_table PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_zetam1_table
  quadmath
  mpfr gmp)

add_executable(build_zeta_trig_tables build_zeta_trig_tables.cpp)
target_include_directories(build_zeta_trig_tables PRIVATE
  $ENV{HOME}/mpreal
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(build_zeta_trig_tables
  quadmath
  mpfr gmp)

add_executable(debug_dirichlet_eta debug_dirichlet_eta.cpp)
target_include_directories(debug_dirichlet_eta PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(debug_dirichlet_eta quadmath)

#add_executable(debug_polylog debug_polylog.cpp)
#target_include_directories(debug_polylog PRIVATE
#  ../../3rdparty)
#target_link_libraries(debug_polylog )

add_executable(test_bose_einstein test_bose_einstein.cpp)
target_include_directories(test_bose_einstein PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_bose_einstein quadmath)

add_executable(test_clausen test_clausen.cpp)
target_include_directories(test_clausen PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_clausen
  quadmath
  -L../../wrappers/debug -lwrap_gsl)

add_executable(test_debye test_debye.cpp)
target_include_directories(test_debye PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_debye
  quadmath
  -L../../wrappers/debug -lwrap_gsl)

add_executable(test_dilog test_dilog.cpp)
target_include_directories(test_dilog PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_dilog quadmath)

add_executable(test_dirichlet_eta test_dirichlet_eta.cpp)
target_include_directories(test_dirichlet_eta PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_dirichlet_eta quadmath)

add_executable(test_fermi_dirac test_fermi_dirac.cpp)
target_include_directories(test_fermi_dirac PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_fermi_dirac
  quadmath
  -L../../wrappers/debug -lwrap_gsl)

add_executable(test_hurwitz_zeta test_hurwitz_zeta.cpp)
target_include_directories(test_hurwitz_zeta PRIVATE
  ../../include
  ../../cxx_math_const/include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_hurwitz_zeta quadmath)

add_executable(test_hurwitz_zeta_new test_hurwitz_zeta_new.cpp)
target_include_directories(test_hurwitz_zeta_new PRIVATE
  ../../include
  ../../cxx_math_const/include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_hurwitz_zeta_new quadmath)

add_executable(test_lerch test_lerch.cpp)
target_include_directories(test_lerch PRIVATE
  $ENV{HOME}/mpreal
  ../.. # For statistics.h
  ../../cxx_summation/include
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_lerch
  mpfr gmp
  quadmath
  -L../../wrappers/debug -lwrap_lerchphi)

add_executable(test_polylog test_polylog.cpp)
target_include_directories(test_polylog PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_polylog
  quadmath
  -L../../wrappers/debug -lwrap_cephes)

add_executable(test_riemann_zeta test_riemann_zeta.cpp)
target_include_directories(test_riemann_zeta PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_riemann_zeta
  quadmath
  -L../../wrappers/debug -lwrap_gsl)

add_executable(zeta_glob_deathmatch zeta_glob_deathmatch.cpp)
target_include_directories(zeta_glob_deathmatch PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(zeta_glob_deathmatch
  quadmath
  -L../../wrappers/debug -lwrap_gsl)
