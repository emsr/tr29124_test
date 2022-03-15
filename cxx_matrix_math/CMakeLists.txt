cmake_minimum_required (VERSION 3.10)

project(
  cxx_matrix_math
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_matrix_math_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_matrix_math INTERFACE)
target_include_directories(cxx_matrix_math INTERFACE include)

add_executable(test_matrix test/src/test_matrix.cpp)
target_link_libraries(test_matrix cxx_matrix_math)
add_test(NAME run_test_matrix COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_matrix > output/test_matrix.txt")

add_executable(test_vandermonde test/src/test_vandermonde.cpp)
target_link_libraries(test_vandermonde cxx_matrix_math cxx_integration)
add_test(NAME run_test_vandermonde COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_vandermonde > output/test_vandermonde.txt")

add_executable(test_cholesky_decomposition test/src/test_cholesky_decomposition.cpp)
target_link_libraries(test_cholesky_decomposition cxx_matrix_math cxx_integration)
add_test(NAME run_test_cholesky_decomposition COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cholesky_decomposition > output/test_cholesky_decomposition.txt")

add_executable(test_lu_decomposition test/src/test_lu_decomposition.cpp)
target_link_libraries(test_lu_decomposition cxx_matrix_math cxx_integration)
add_test(NAME run_test_lu_decomposition COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lu_decomposition > output/test_lu_decomposition.txt")

add_executable(test_sv_decomposition test/src/test_sv_decomposition.cpp)
target_link_libraries(test_sv_decomposition cxx_matrix_math cxx_integration)
add_test(NAME run_test_sv_decomposition COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sv_decomposition > output/test_sv_decomposition.txt")
