cmake_minimum_required (VERSION 3.10)

project(
  gst_scorer
  VERSION 1.0
  LANGUAGES Fortran
)

set(src_dir Fortran77/Src/Dp)
set(test_dir Fortran77/Drivers/Dp)

add_library(gst_scorer SHARED ${src_dir}/src_819.f ${src_dir}/scorer.f90)
target_link_libraries(gst_scorer gst_airy)

# I think the file names were switched by the authors: src_819.f <-> src.f
# This src.f is exactly the same as 819_Airy/Fortran77/Src/Dp/src.f
# ${src_dir}/src.f

add_executable(sctest ${test_dir}/driver.f)
target_link_libraries(sctest gst_scorer)

add_executable(test_gst_scorer ${src_dir}/test_scorer.cpp)
target_link_libraries(test_gst_scorer gst_scorer)
