cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_summation
  VERSION 0.5
  LANGUAGES CXX Fortran
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

add_custom_target(make_lab_summation_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_cohen_summation test_cohen_summation.cpp)
target_link_libraries(test_cohen_summation ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_cohen_summation COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cohen_summation > output/test_cohen_summation.txt")

add_executable(test_CohenSum test_CohenSum.cpp)
target_link_libraries(test_CohenSum ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_CohenSum COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_CohenSum > output/test_CohenSum.txt")
