cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_norm_functions
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

add_custom_target(make_lab_norm_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_abs test_abs.cpp)
add_test(NAME run_test_abs COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_abs > output/test_abs.txt")

add_executable(test_arith_geom_mean test_arith_geom_mean.cpp)
target_link_libraries(test_arith_geom_mean
  ${TR29124_LIBRARIES}
  cxx_float128
  quadmath
  wrap_boost)
add_test(NAME run_test_arith_geom_mean COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_arith_geom_mean > output/test_arith_geom_mean.txt")

add_executable(test_hypot test_hypot.cpp)
target_link_libraries(test_hypot
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_hypot COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hypot > output/test_hypot.txt")

add_executable(test_logsumexp test_logsumexp.cpp)
add_test(NAME run_test_logsumexp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_logsumexp > output/test_logsumexp.txt")

add_executable(test_power_mean test_power_mean.cpp)
add_test(NAME run_test_power_mean COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_power_mean > output/test_power_mean.txt")

add_executable(test_power_norm test_power_norm.cpp)
add_test(NAME run_test_power_norm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_power_norm > output/test_power_norm.txt")
