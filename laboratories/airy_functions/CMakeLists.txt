cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")

enable_testing()

add_executable(airy_toy airy_toy.cpp)
target_include_directories(airy_toy PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(airy_toy quadmath)

add_executable(airy_toy_old airy_toy_old.cpp)
target_include_directories(airy_toy_old PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(airy_toy_old quadmath)

add_executable(debug_airy debug_airy.cpp)
target_include_directories(debug_airy PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(debug_airy quadmath -L../../wrappers/debug -lwrap_gsl)

add_executable(plot_airy plot_airy.cpp)
target_include_directories(plot_airy PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(plot_airy quadmath)

add_executable(test_airy_roots test_airy_roots.cpp)
target_include_directories(test_airy_roots PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(test_airy_roots quadmath)

add_executable(test_little_airy test_little_airy.cpp)
target_include_directories(test_little_airy PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_little_airy quadmath -L../../wrappers/debug -lwrap_gsl)
