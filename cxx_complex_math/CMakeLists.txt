cmake_minimum_required (VERSION 3.10)

project(
  cxx_complex_math
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")

enable_testing()

add_custom_target(make_cxx_complex_math_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_complex_math INTERFACE)
target_include_directories(cxx_complex_math INTERFACE include)
target_link_libraries(cxx_complex_math INTERFACE cxx_numeric_limits cxx_math_constants)

add_executable(test_complex_safe_math test_complex_safe_math.cpp)
target_link_libraries(test_complex_safe_math cxx_complex_math cxx_math_constants)

add_executable(test_complex_norms test_complex_norms.cpp)
target_link_libraries(test_complex_norms cxx_complex_math)

add_executable(test_complex_midpoint test_complex_midpoint.cpp)
target_link_libraries(test_complex_midpoint cxx_complex_math)

add_executable(test_complex_lerp test_complex_lerp.cpp)
target_link_libraries(test_complex_lerp cxx_complex_math)
