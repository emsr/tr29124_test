cmake_minimum_required (VERSION 3.10)

project(
  gst_airy
  VERSION 1.0
  LANGUAGES Fortran
)

set(lib_src_dir Fortran77/Dp/Src)
set(test_src_dir Fortran77/Dp/Drivers)

add_library(gst_airy SHARED ${lib_src_dir}/port.f ${lib_src_dir}/src.f)

add_executable(test_gst_airy ${test_src_dir}/driver.f)
target_link_libraries(test_gst_airy gst_airy)
