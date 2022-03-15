cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_floating_point_tools
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

add_custom_target(make_lab_floating_point_tools_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(builtins builtins.cpp)
add_test(NAME run_builtins COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/builtins > output/builtins.txt")

add_executable(test_float128 test_float128.cpp)
target_link_libraries(test_float128 ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_float128 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_float128 > output/test_float128.txt")

add_executable(test_limits test_limits.cpp)
target_link_libraries(test_limits ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_limits COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_limits > output/test_limits.txt")

if (HAVE_MPREAL)
  add_executable(test_maxint test_maxint.cpp)
  target_link_libraries(test_maxint ${EMSR_LIBRARIES} quadmath mpreal gmp mpfr)
  add_test(NAME run_test_maxint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_maxint > output/test_maxint.txt")
endif(HAVE_MPREAL)

add_executable(test_mod2pi test_mod2pi.cpp)
target_link_libraries(test_mod2pi ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_mod2pi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_mod2pi > output/test_mod2pi.txt")

add_executable(test_lab_ulp test_ulp.cpp)
add_test(NAME run_test_lab_ulp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lab_ulp > output/test_lab_ulp.txt")

add_executable(variant_help variant_help.cpp)
add_test(NAME run_variant_help COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/variant_help > output/variant_help.txt")
