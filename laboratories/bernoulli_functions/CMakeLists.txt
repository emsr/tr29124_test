cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_bernoulli_functions
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

add_custom_target(make_lab_bernoulli_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(build_bernoulli_2n_table build_bernoulli_2n_table.cpp)
target_link_libraries(build_bernoulli_2n_table ${TR29124_LIBRARIES} quadmath mpreal mpfr)
add_test(NAME run_build_bernoulli_2n_table COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/build_bernoulli_2n_table > output/build_bernoulli_2n_table.txt")

add_executable(test_bernoulli test_bernoulli.cpp)
target_include_directories(test_bernoulli PRIVATE
  ../../cxx_float128/include)
target_link_libraries(test_bernoulli
  ${TR29124_LIBRARIES}
  cxx_float128
  quadmath
  wrap_burkhardt
  gfortran)
add_test(NAME run_test_bernoulli COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bernoulli > output/test_bernoulli.txt")
