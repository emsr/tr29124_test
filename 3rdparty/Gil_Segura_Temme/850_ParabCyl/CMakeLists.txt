cmake_minimum_required (VERSION 3.10)

project(
  gst_parab_cyl
  VERSION 1.0
  LANGUAGES Fortran
)

set(lib_src_dir Fortran90/Dp/Src)
set(lib_src2_dir Fortran90/Dp/Src2.0)
set(test_src_dir Fortran90/Dp/Drivers)

add_library(gst_parab_cyl SHARED ${lib_src_dir}/air.f90 ${lib_src_dir}/cons.f90 ${lib_src_dir}/parab.f90)
add_library(gst_parab_cyl2 SHARED ${lib_src2_dir}/air.f90 ${lib_src_dir}/cons.f90 ${lib_src_dir}/parab.f90)

add_executable(test_gst_parab_cyl ${test_src_dir}/driver.f90)
target_link_libraries(test_gst_parab_cyl gst_parab_cyl)

add_executable(test_gst_parab_cyl2 ${test_src_dir}/driver.f90)
target_link_libraries(test_gst_parab_cyl2 gst_parab_cyl2)
