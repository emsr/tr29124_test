cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_basic_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_lab_basic_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_q_number test_q_number.cpp)
add_test(NAME run_test_q_number COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_q_number > output/test_q_number.txt")

add_executable(test_q_pochhammer test_q_pochhammer.cpp)
add_test(NAME run_test_q_pochhammer COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_q_pochhammer > output/test_q_pochhammer.txt")
