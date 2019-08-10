cmake_minimum_required (VERSION 3.10)

project(
  toms917_wright_omega
  VERSION 1.0
  LANGUAGES C
)

enable_testing()

set(lib_src_dir C/Src)
set(test_src_dir C/Drivers)

add_library(toms917_wright_omega SHARED ${lib_src_dir}/wright.c)
target_link_libraries(toms917_wright_omega m)
target_include_directories(toms917_wright_omega PUBLIC ${lib_src_dir})

add_executable(test_toms917_wright_omega ${test_src_dir}/driver.c)
target_link_libraries(test_toms917_wright_omega toms917_wright_omega)

add_executable(test_toms917_wright_omega_boundary ${test_src_dir}/test_boundary.c)
target_link_libraries(test_toms917_wright_omega_boundary toms917_wright_omega)
