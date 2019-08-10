cmake_minimum_required (VERSION 3.10)

project(
  toms804_mathieu
  VERSION 1.0
  LANGUAGES CXX
)

enable_testing()

set(lib_src_dir CPP/Dp/Src)
set(test_src_dir CPP/Dp/Drivers)

add_library(toms804_mathieu SHARED
  ${lib_src_dir}/bsslr.cpp
  ${lib_src_dir}/mathur.cpp
  ${lib_src_dir}/mbsslr.cpp
  ${lib_src_dir}/mcnr.cpp
  ${lib_src_dir}/mmathur.cpp)
target_include_directories(toms804_mathieu PUBLIC ${lib_src_dir})

add_executable(test_toms804_mathieu ${test_src_dir}/driver.cpp)
target_link_libraries(test_toms804_mathieu toms804_mathieu)
