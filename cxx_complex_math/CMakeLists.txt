cmake_minimum_required (VERSION 3.10)

project(
  cxx_complex_math
  VERSION 0.5
  LANGUAGES CXX
)

enable_testing()

add_library(cxx_complex_math INTERFACE)
target_include_directories(cxx_complex_math INTERFACE include)
