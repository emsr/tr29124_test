cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_coulomb_functions
  VERSION 0.5
  LANGUAGES CXX Fortran
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_executable(run_coulfg run_coulfg.cpp coulfg.cpp)
add_test(NAME run_run_coulfg COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/run_coulfg > run_coulfg.txt")

add_executable(test_coulomb test_coulomb.cpp)
add_test(NAME run_test_coulomb COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_coulomb > test_coulomb.txt")

add_executable(test_hydrogen test_hydrogen.cpp)
target_include_directories(test_hydrogen PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_hydrogen quadmath)
add_test(NAME run_test_hydrogen COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hydrogen > test_hydrogen.txt")
