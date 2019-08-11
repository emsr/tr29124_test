cmake_minimum_required (VERSION 3.10)

project(
  gst_parabolic_cylinder
  VERSION 1.0
  LANGUAGES Fortran
)

set(lib_src_dir Fortran90/Dp/Src)
set(lib_src2_dir Fortran90/Dp/Src2.0)
set(test_src_dir Fortran90/Dp/Drivers)

add_library(gst_parabolic_cylinder SHARED
  ${lib_src_dir}/air.f90
  ${lib_src_dir}/cons.f90
  ${lib_src_dir}/parab.f90)
#target_include_directories(gst_parabolic_cylinder PUBLIC .)

add_executable(test_gst_parab_cyl1 ${test_src_dir}/driver.f90)
target_link_libraries(test_gst_parab_cyl1 gst_parabolic_cylinder)

# Version 2 is the real goods.

add_library(gst_parabolic_cylinder_2 SHARED
  ${lib_src2_dir}/air.f90
  ${lib_src2_dir}/cons.f90
  ${lib_src2_dir}/parab.f90
  ${lib_src2_dir}/parabolic_cylinder.f90)
target_include_directories(gst_parabolic_cylinder_2 PUBLIC .)

add_executable(test_gst_parab_cyl2 ${test_src_dir}/driver.f90)
target_link_libraries(test_gst_parab_cyl2 gst_parabolic_cylinder_2)

add_executable(test_gst_parabolic_cylinder ${lib_src2_dir}/test_parabolic_cylinder.cpp)
target_link_libraries(test_gst_parabolic_cylinder gst_parabolic_cylinder_2)
