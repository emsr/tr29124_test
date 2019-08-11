cmake_minimum_required (VERSION 3.10)

project(
  gst_incomplete_gamma
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_incomplete_gamma SHARED cons.f90 incgam.f90 incomplete_gamma.f90)
target_include_directories(gst_incomplete_gamma PUBLIC .)

add_executable(test_incgam testincgam.f90)
target_link_libraries(test_incgam gst_incomplete_gamma)

add_executable(test_gst_incomplete_gamma test_inc_gamma.cpp)
target_link_libraries(test_gst_incomplete_gamma gst_incomplete_gamma)
