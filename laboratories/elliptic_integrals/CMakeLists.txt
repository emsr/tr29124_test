cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_elliptic_integrals
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

#add_executable(carlson1 carlson1.cpp)
#target_include_directories(carlson1 PRIVATE )
#target_link_libraries(carlson1 )

#add_executable(carlson2 carlson2.cpp)
#target_include_directories(carlson2 PRIVATE )
#target_link_libraries(carlson2 )

#add_executable(carlson3 carlson3.cpp)
#target_include_directories(carlson3 PRIVATE )
#target_link_libraries(carlson3 )

#add_executable(carlson4 carlson4.cpp)
#target_include_directories(carlson4 PRIVATE )
#target_link_libraries(carlson4 )

#add_executable(carlson carlson.cpp)
#target_include_directories(carlson PRIVATE )
#target_link_libraries(carlson )

add_executable(carlson_std carlson_std.cpp)
target_include_directories(carlson_std PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(carlson_std quadmath)
add_test(NAME run_carlson_std COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/carlson_std > carlson_std.txt")

add_executable(test_comp_ellint_1 test_comp_ellint_1.cpp)
target_include_directories(test_comp_ellint_1 PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_comp_ellint_1 quadmath)
add_test(NAME run_test_comp_ellint_1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_comp_ellint_1 > test_comp_ellint_1.txt")

add_executable(test_heuman_lambda test_heuman_lambda.cpp)
target_include_directories(test_heuman_lambda PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_heuman_lambda
  quadmath
  -L../../wrappers -lwrap_boost)
add_test(NAME run_test_heuman_lambda COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_heuman_lambda > test_heuman_lambda.txt")

add_executable(test_jacobi_zeta test_jacobi_zeta.cpp)
target_include_directories(test_jacobi_zeta PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_jacobi_zeta quadmath -L../../wrappers -lwrap_boost)
add_test(NAME run_test_jacobi_zeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_zeta > test_jacobi_zeta.txt")

add_executable(test_legendre_ellint test_legendre_ellint.cpp)
target_include_directories(test_legendre_ellint PRIVATE
  ../../wrappers
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_legendre_ellint
  quadmath
  -L../../wrappers -lwrap_gsl)
add_test(NAME run_test_legendre_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre_ellint > test_legendre_ellint.txt")
