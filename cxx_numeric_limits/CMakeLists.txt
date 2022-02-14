cmake_minimum_required (VERSION 3.10)

project(
  cxx_numeric_limits
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_numeric_limits_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_numeric_limits INTERFACE)
target_include_directories(cxx_numeric_limits INTERFACE include)
target_link_libraries(cxx_numeric_limits INTERFACE mpreal multiprecision)

if (HAVE_MPREAL)
  add_executable(test_numeric_limits test_numeric_limits.cpp)
  target_link_libraries(test_numeric_limits cxx_numeric_limits cxx_fp_utils cxx_float128 quadmath mpreal gmp mpfr)
  add_test(NAME run_test_numeric_limits COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_numeric_limits > output/test_numeric_limits.txt")
endif(HAVE_MPREAL)
