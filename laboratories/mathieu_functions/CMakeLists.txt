cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_mathieu_functions
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

add_executable(test_mathieu test_mathieu.cpp)
target_include_directories(test_mathieu PRIVATE )
target_link_libraries(test_mathieu quadmath)
add_test(NAME run_test_mathieu COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_mathieu > output/test_mathieu.txt")
