cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_theta_functions
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

add_executable(test_jacobi_ellint test_jacobi_ellint.cpp)
target_include_directories(test_jacobi_ellint PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_jacobi_ellint quadmath -L../../wrappers -lwrap_gsl)
add_test(NAME run_test_jacobi_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_ellint > output/test_jacobi_ellint.txt")

add_executable(test_jacobi_inv test_jacobi_inv.cpp)
target_include_directories(test_jacobi_inv PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_jacobi_inv quadmath)
add_test(NAME run_test_jacobi_inv COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_inv > output/test_jacobi_inv.txt")

add_executable(test_jacobi_theta test_jacobi_theta.cpp)
target_include_directories(test_jacobi_theta PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_jacobi_theta quadmath -L../../wrappers -lwrap_gsl)
add_test(NAME run_test_jacobi_theta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_theta > output/test_jacobi_theta.txt")

add_executable(test_theta test_theta.cpp)
target_include_directories(test_theta PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_theta quadmath)
add_test(NAME run_test_theta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_theta > output/test_theta.txt")

add_executable(test_weierstrass_ellint test_weierstrass_ellint.cpp)
target_include_directories(test_weierstrass_ellint PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_weierstrass_ellint quadmath)
add_test(NAME run_test_weierstrass_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_weierstrass_ellint > output/test_weierstrass_ellint.txt")
