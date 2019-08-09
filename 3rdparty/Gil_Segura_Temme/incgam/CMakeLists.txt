cmake_minimum_required (VERSION 3.10)

project(
  gst_incgam
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_incgam SHARED cons.f90 incgam.f90 incomplete_gamma.f90)

add_executable(test_incgam testincgam.f90)
target_link_libraries(test_incgam gst_incgam)

add_executable(test_gst_inc_gamma test_inc_gamma.cpp)
target_link_libraries(test_gst_inc_gamma gst_incgam)
