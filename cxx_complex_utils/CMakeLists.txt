cmake_minimum_required (VERSION 3.10)

project(
  cxx_complex_utils
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_complex_utils_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_complex_utils INTERFACE)
target_include_directories(cxx_complex_utils INTERFACE include)

add_executable(test_is_complex test_is_complex.cpp)
target_include_directories(test_is_complex PRIVATE
  ../include)
target_link_libraries(test_is_complex cxx_complex_utils cxx_fp_utils quadmath)
