cmake_minimum_required (VERSION 3.10)

project(
  cxx_maths
  VERSION 0.5
  LANGUAGES CXX
)

enable_testing()

add_subdirectory(wrappers)
add_subdirectory(libstdc++_support)
add_subdirectory(laboratory)
