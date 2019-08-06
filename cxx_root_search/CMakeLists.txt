cmake_minimum_required (VERSION 3.10)

project(
  cxx_root_search
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_cxx_root_search_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_root_search INTERFACE)
target_include_directories(cxx_root_search
  INTERFACE include)

add_executable(test_root_search test_root_search.cpp)
target_include_directories(test_root_search PRIVATE
  .)
target_link_libraries(test_root_search
  cxx_math_constants
  cxx_polynomial
  cxx_root_search)
add_test(NAME test_root_search COMMAND bash -c "./bin/test_root_search > test_root_search.txt")
