cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_theta_functions
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

add_custom_target(make_lab_theta_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_jacobi_ellint test_jacobi_ellint.cpp)
target_link_libraries(test_jacobi_ellint ${EMSR_LIBRARIES} quadmath wrap_gsl)
add_test(NAME run_test_jacobi_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_ellint > output/test_jacobi_ellint.txt")

add_executable(test_jacobi_inv test_jacobi_inv.cpp)
target_link_libraries(test_jacobi_inv ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_jacobi_inv COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_inv > output/test_jacobi_inv.txt")

add_executable(test_jacobi_theta test_jacobi_theta.cpp)
target_include_directories(test_jacobi_theta PRIVATE
  ../../wrappers)
target_link_libraries(test_jacobi_theta ${EMSR_LIBRARIES} cxx_float128 quadmath wrap_gsl)
add_test(NAME run_test_jacobi_theta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_theta > output/test_jacobi_theta.txt")

add_executable(test_theta test_theta.cpp)
target_link_libraries(test_theta ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_theta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_theta > output/test_theta.txt")

add_executable(test_weierstrass_ellint test_weierstrass_ellint.cpp)
target_link_libraries(test_weierstrass_ellint ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_weierstrass_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_weierstrass_ellint > output/test_weierstrass_ellint.txt")
