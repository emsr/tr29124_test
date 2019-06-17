cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_recursive_functions
  VERSION 0.1
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

add_executable(test_fibonacci test_fibonacci.cpp)
target_include_directories(test_fibonacci PRIVATE
  ../../cxx_math_constants/include
  ../../cxx_fp_utils/include
  ../../include # All this rest for cos_pi.
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_fibonacci quadmath)
add_test(NAME run_test_fibonacci COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_fibonacci > test_fibonacci.txt")
