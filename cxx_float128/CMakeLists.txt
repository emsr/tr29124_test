cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_float128
  VERSION 0.1
  LANGUAGES C CXX
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

# Header-only library.
add_library(cxx_float128 INTERFACE)
target_include_directories(cxx_float128 INTERFACE include)
