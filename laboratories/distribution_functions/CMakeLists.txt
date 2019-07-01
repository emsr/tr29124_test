cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_distribution_functions
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

add_executable(test_marcum_q test_marcum_q.cpp)
target_include_directories(test_marcum_q PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../polynomial/include
  ../../cxx_continued_fractions/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_marcum_q quadmath)
add_test(NAME run_test_marcum_q COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_marcum_q > output/test_marcum_q.txt")
