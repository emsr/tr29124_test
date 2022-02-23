cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_polunomials
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_lab_polynomials_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(rational_fit rational_fit.cpp)
target_link_libraries(rational_fit cxx_matrix_math cxx_polynomial)
add_test(NAME rational_fit COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/rational_fit > output/rational_fit.txt")

add_executable(rational_fit_new rational_fit_new.cpp)
target_link_libraries(rational_fit_new cxx_matrix_math cxx_polynomial)
add_test(NAME rational_fit_new COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/rational_fit_new > output/rational_fit_new.txt")

add_executable(rational_fit_old rational_fit_old.cpp)
target_link_libraries(rational_fit_old cxx_matrix_math cxx_polynomial)
add_test(NAME rational_fit_old COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/rational_fit_old > output/rational_fit_old.txt")

add_executable(rational_fit_help help.cpp)
target_link_libraries(rational_fit_help cxx_complex_utils cxx_integration cxx_special_functions)
