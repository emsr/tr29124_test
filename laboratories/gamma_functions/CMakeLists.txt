cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_gamma_functions
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

add_executable(build_gamma_lanczos build_gamma_lanczos.cpp)
target_include_directories(build_gamma_lanczos PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(build_gamma_lanczos quadmath)

add_executable(build_gamma_recip build_gamma_recip.cpp)
target_include_directories(build_gamma_recip PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  $ENV{HOME}/mpreal
)
target_link_libraries(build_gamma_recip quadmath mpfr gmp)

add_executable(build_gamma_spouge build_gamma_spouge.cpp)
target_include_directories(build_gamma_spouge PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(build_gamma_spouge quadmath)

add_executable(plot_gamma plot_gamma.cpp)
target_include_directories(plot_gamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(plot_gamma quadmath)

add_executable(test_binet test_binet.cpp)
target_include_directories(test_binet PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_rational/include
)
target_link_libraries(test_binet quadmath)

add_executable(test_binet_float test_binet_float.cpp)
target_include_directories(test_binet_float PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(test_binet_float quadmath)

add_executable(test_complex_gamma test_complex_gamma.cpp)
target_include_directories(test_complex_gamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(test_complex_gamma quadmath)

add_executable(test_digamma test_digamma.cpp)
target_include_directories(test_digamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_digamma quadmath -L../../wrappers/debug -lwrap_gsl)

add_executable(test_factorial test_factorial.cpp)
target_include_directories(test_factorial PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(test_factorial quadmath)

add_executable(test_falling_factorial test_falling_factorial.cpp)
target_include_directories(test_falling_factorial PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_falling_factorial quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_gamma test_gamma.cpp)
target_include_directories(test_gamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_gamma quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_gamma_ratio test_gamma_ratio.cpp)
target_include_directories(test_gamma_ratio PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_gamma_ratio quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_gamma_reciprocal test_gamma_reciprocal.cpp)
target_include_directories(test_gamma_reciprocal PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_gamma_reciprocal quadmath)

add_executable(test_inv_gamma test_inv_gamma.cpp)
target_include_directories(test_inv_gamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_inv_gamma quadmath)

add_executable(test_inv_lgamma test_inv_lgamma.cpp)
target_include_directories(test_inv_lgamma PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_inv_lgamma quadmath)

add_executable(test_polygamma test_polygamma.cpp)
target_include_directories(test_polygamma PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_continued_fractions/include
  ../../wrappers)
target_link_libraries(test_polygamma quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_rising_factorial test_rising_factorial.cpp)
target_include_directories(test_rising_factorial PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_rising_factorial quadmath -L../../wrappers/debug -lwrap_boost)
