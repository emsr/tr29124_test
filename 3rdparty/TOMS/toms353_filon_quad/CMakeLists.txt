cmake_minimum_required (VERSION 3.10)

project(
  toms353_filon_quad
  VERSION 1.0
  LANGUAGES Fortran
)

enable_testing()

set(lib_src_dir Fortran/Sp/Src)
set(test_src_dir Fortran/Sp/Drivers)

add_library(toms353_filon_quad SHARED ${lib_src_dir}/src.f)

add_executable(test_toms353_filon_quad ${test_src_dir}/driver.f)
target_link_libraries(test_toms353_filon_quad toms353_filon_quad)
