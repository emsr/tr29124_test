cmake_minimum_required (VERSION 3.10)

project(
  toms804_adapt_dubdiv_quad
  VERSION 1.0
  LANGUAGES Fortran
)

enable_testing()

set(lib_src_dir Fortran/Sp/Src)
set(test_src_dir Fortran/Sp/Drivers)

add_library(toms804_adapt_dubdiv_quad SHARED ${lib_src_dir}/src.f)

add_executable(test_toms804_adapt_dubdiv_quad ${test_src_dir}/driver.f)
target_link_libraries(test_toms804_adapt_dubdiv_quad toms804_adapt_dubdiv_quad)
