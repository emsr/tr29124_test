cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_float128
  VERSION 0.1
  LANGUAGES C CXX
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

add_custom_target(make_cxx_float128_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

# Header-only library.
add_library(cxx_float128 INTERFACE)
target_link_libraries(cxx_float128 INTERFACE cxx_numeric_limits)
target_include_directories(cxx_float128 INTERFACE include)
