cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_beta_functions
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

add_custom_target(make_lab_beta_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_beta test_beta.cpp)
target_link_libraries(test_beta ${EMSR_LIBRARIES} cxx_float128
  quadmath
  wrap_boost)
add_test(NAME run_test_beta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_beta > output/test_beta.txt")

add_executable(test_beta_inc test_beta_inc.cpp)
target_link_libraries(test_beta_inc ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_beta_inc COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_beta_inc > output/test_beta_inc.txt")

add_executable(test_inv_ibeta test_inv_ibeta.cpp)
target_include_directories(test_inv_ibeta PRIVATE
  ../../cxx_root_search/include)
target_link_libraries(test_inv_ibeta ${EMSR_LIBRARIES} quadmath)
# FIXME add_test(NAME run_test_inv_ibeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_inv_ibeta > output/test_inv_ibeta.txt")
