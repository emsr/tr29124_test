cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_interval
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_interval_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_interval test_interval.cpp)
target_include_directories(test_interval PRIVATE include)
add_test(NAME run_test_interval COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_interval > output/test_interval.txt")
