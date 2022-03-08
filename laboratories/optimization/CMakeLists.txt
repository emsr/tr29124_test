cmake_minimum_required (VERSION 3.18)

include(CheckCXXCompilerFlag)

project(
  lab_optimization
  VERSION 0.2
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(lab_optimization INTERFACE)
target_include_directories(lab_optimization INTERFACE .)

add_executable(test_lab_optimization test_optimization.cpp)
target_link_libraries(test_lab_optimization lab_optimization)
