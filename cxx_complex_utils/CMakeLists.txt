cmake_minimum_required (VERSION 3.10)

project(
  cxx_complex_utils
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_executable(test_is_complex test_is_complex.cpp)
target_include_directories(test_is_complex PRIVATE
  ../cxx_fp_utils/include
  include
  ../include)
target_link_libraries(test_is_complex quadmath)
