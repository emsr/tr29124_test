cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_special_functions
  VERSION 0.9
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_special_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_special_functions INTERFACE)
target_include_directories(cxx_special_functions INTERFACE ../include)
target_link_libraries(cxx_special_functions INTERFACE cxx_complex_math cxx_numeric_limits)
