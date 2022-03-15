cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_parabolic_cylinder_functions
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

add_custom_target(make_lab_parabolic_cylinder_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_parab_cyl test_parab_cyl.cpp)
target_include_directories(test_parab_cyl PRIVATE .)
target_link_libraries(test_parab_cyl ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_parab_cyl COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_parab_cyl > output/test_parab_cyl.txt")
