cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  multiprecision
  VERSION 0.1
  LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_library(multiprecision INTERFACE)
target_include_directories(multiprecision INTERFACE .)

if (HAVE_MPREAL)
  add_executable(mpfr_gexpr mpfr_gexpr.c)
  target_link_libraries(mpfr_gexpr mpfr gmp m)
else()
  message(WARNING "No MPReal!")
endif(HAVE_MPREAL)
