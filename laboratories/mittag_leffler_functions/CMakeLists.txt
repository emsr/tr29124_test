cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_mittag_leffler_functions
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

add_custom_target(make_lab_mittag_leffler_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_mittag_leffler test_mittag_leffler.cpp)
target_link_libraries(test_mittag_leffler ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_mittag_leffler COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_mittag_leffler > output/test_mittag_leffler.txt")
