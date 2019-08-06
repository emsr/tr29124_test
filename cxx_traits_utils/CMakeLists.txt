cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_traits_utils
  VERSION 0.1.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_cxx_traits_utils_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_traits_utils INTERFACE)
target_include_directories(cxx_traits_utils INTERFACE include)

add_executable(test_traits_utils test_traits_utils.cpp)
target_link_libraries(test_traits_utils PRIVATE cxx_traits_utils)

add_test(NAME run_test_traits_utils COMMAND bash -c "./bin/test_traits_utils > output/test_traits_utils.txt")
