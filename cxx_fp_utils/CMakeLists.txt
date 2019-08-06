cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_fp_utils
  VERSION 1.0.0
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

add_custom_target(make_cxx_fp_utils_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_fp_utils INTERFACE)
target_include_directories(cxx_fp_utils INTERFACE include)

add_executable(test_fp_utils test_fp_utils.cpp)
target_link_libraries(test_fp_utils cxx_fp_utils)
add_test(NAME run_test_fp_utils COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_fp_utils > output/test_fp_utils.txt")

add_executable(test_ulp test_ulp.cpp)
target_link_libraries(test_ulp cxx_fp_utils)
add_test(NAME run_test_ulp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_ulp > output/test_ulp.txt")
