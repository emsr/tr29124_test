cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_differentiation
  VERSION 1.0.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

list(APPEND cxx_differentiation_INCLUDES
  include/ext/derivative.h
  include/ext/derivative.tcc
)


add_executable(test_derivative test_derivative.cpp)
target_include_directories(test_derivative PRIVATE include)

add_test(NAME derivative_test0 COMMAND bash -c "./test_derivative > test_derivative.txt")
