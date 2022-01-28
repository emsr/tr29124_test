cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_zeta_functions
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

add_custom_target(make_lab_zeta_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(build_etam1_table build_etam1_table.cpp)
target_link_libraries(build_etam1_table
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_etam1_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_etam1_table > output/build_etam1_table.txt")

add_executable(build_nfact_zetanp1 build_nfact_zetanp1.cpp)
target_link_libraries(build_nfact_zetanp1
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_nfact_zetanp1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_nfact_zetanp1 > output/build_nfact_zetanp1.txt")

add_executable(build_zeta_deriv_table build_zeta_deriv_table.cpp)
target_link_libraries(build_zeta_deriv_table
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_zeta_deriv_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_zeta_deriv_table > output/build_zeta_deriv_table.txt")

add_executable(build_zetahalfm1_table build_zetahalfm1_table.cpp)
target_link_libraries(build_zetahalfm1_table
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_zetahalfm1_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_zetahalfm1_table > output/build_zetahalfm1_table.txt")

add_executable(build_zetam1_table build_zetam1_table.cpp)
target_link_libraries(build_zetam1_table
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_zetam1_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_zetam1_table > output/build_zetam1_table.txt")

add_executable(build_zeta_trig_tables build_zeta_trig_tables.cpp)
target_link_libraries(build_zeta_trig_tables
  ${EMSR_LIBRARIES}
  quadmath
  mpreal mpfr gmp)
add_test(NAME run_build_zeta_trig_tables COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_zeta_trig_tables > output/build_zeta_trig_tables.txt")

add_executable(debug_dirichlet_eta debug_dirichlet_eta.cpp)
target_link_libraries(debug_dirichlet_eta ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_debug_dirichlet_eta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/debug_dirichlet_eta > output/debug_dirichlet_eta.txt")

#add_executable(debug_polylog debug_polylog.cpp)
#target_include_directories(debug_polylog PRIVATE
#  ../../3rdparty)
#target_link_libraries(debug_polylog ${EMSR_LIBRARIES} )
#add_test(NAME run_debug_polylog COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/debug_polylog > output/debug_polylog.txt")

add_executable(test_bose_einstein test_bose_einstein.cpp)
target_include_directories(test_bose_einstein PRIVATE
  ../../cxx_float128/include)
target_link_libraries(test_bose_einstein ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_bose_einstein COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bose_einstein > output/test_bose_einstein.txt")

add_executable(test_clausen test_clausen.cpp)
target_link_libraries(test_clausen
  ${EMSR_LIBRARIES}
  quadmath
  wrap_gsl)
add_test(NAME run_test_clausen COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_clausen > output/test_clausen.txt")

add_executable(test_debye test_debye.cpp)
target_link_libraries(test_debye
  ${EMSR_LIBRARIES}
  quadmath
  wrap_gsl)
add_test(NAME run_test_debye COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_debye > output/test_debye.txt")

add_executable(test_dilog test_dilog.cpp)
target_include_directories(test_dilog PRIVATE)
target_link_libraries(test_dilog ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_dilog COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dilog > output/test_dilog.txt")

add_executable(test_dirichlet_eta test_dirichlet_eta.cpp)
target_include_directories(test_dirichlet_eta PRIVATE)
target_link_libraries(test_dirichlet_eta ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_dirichlet_eta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dirichlet_eta > output/test_dirichlet_eta.txt")

add_executable(test_fermi_dirac test_fermi_dirac.cpp)
target_link_libraries(test_fermi_dirac
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  wrap_gsl)
add_test(NAME run_test_fermi_dirac COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_fermi_dirac > output/test_fermi_dirac.txt")

add_executable(test_hurwitz_zeta test_hurwitz_zeta.cpp)
target_link_libraries(test_hurwitz_zeta ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_hurwitz_zeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hurwitz_zeta > output/test_hurwitz_zeta.txt")

add_executable(test_hurwitz_zeta_new test_hurwitz_zeta_new.cpp)
target_link_libraries(test_hurwitz_zeta_new ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_hurwitz_zeta_new COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hurwitz_zeta_new > output/test_hurwitz_zeta_new.txt")

add_executable(test_lerch test_lerch.cpp)
target_include_directories(test_lerch PRIVATE
  ../.. # For statistics.h
)
target_link_libraries(test_lerch
  ${EMSR_LIBRARIES}
  mpreal mpfr gmp
  quadmath
  wrap_lerchphi)
add_test(NAME run_test_lerch COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lerch > output/test_lerch.txt")

add_executable(test_polylog test_polylog.cpp)
target_link_libraries(test_polylog
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  wrap_cephes)
add_test(NAME run_test_polylog COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_polylog > output/test_polylog.txt")

add_executable(test_riemann_zeta test_riemann_zeta.cpp)
target_link_libraries(test_riemann_zeta
  ${EMSR_LIBRARIES}
  cxx_float128
  quadmath
  wrap_gsl)
add_test(NAME run_test_riemann_zeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_riemann_zeta > output/test_riemann_zeta.txt")

add_executable(zeta_glob_deathmatch zeta_glob_deathmatch.cpp)
target_link_libraries(zeta_glob_deathmatch
  cxx_math_constants
  quadmath
  wrap_gsl)
# FIXME: This passes but takes forever and times out on most CI. Put it in another category.
#add_test(NAME run_zeta_glob_deathmatch COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/zeta_glob_deathmatch > output/zeta_glob_deathmatch.txt")
