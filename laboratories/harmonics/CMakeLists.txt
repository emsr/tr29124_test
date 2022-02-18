cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_harmonics
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

add_custom_target(make_lab_harmonics_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(oblate_harmonic oblate_harmonic.cpp)
add_test(NAME run_oblate_harmonic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/oblate_harmonic > output/oblate_harmonic.txt")

add_executable(prolate_harmonic prolate_harmonic.cpp)
add_test(NAME run_prolate_harmonic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/prolate_harmonic > output/prolate_harmonic.txt")

add_executable(toroidal_harmonic toroidal_harmonic.cpp)
add_test(NAME run_toroidal_harmonic COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/toroidal_harmonic > output/toroidal_harmonic.txt")
