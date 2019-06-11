cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_mittag_leffler_functions
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

add_executable(test_mittag_leffler test_mittag_leffler.cpp)
target_include_directories(test_mittag_leffler PRIVATE
  ../../cxx_math_const/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include)
target_link_libraries(test_mittag_leffler quadmath)
add_test(NAME run_test_mittag_leffler COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_mittag_leffler > test_mittag_leffler.txt")
