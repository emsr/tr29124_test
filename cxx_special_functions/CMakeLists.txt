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

add_library(cxx_special_functions INTERFACE)
target_include_directories(cxx_special_functions INTERFACE ../include)
