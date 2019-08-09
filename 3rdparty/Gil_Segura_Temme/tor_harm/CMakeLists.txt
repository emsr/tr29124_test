cmake_minimum_required (VERSION 3.10)

project(
  gst_tor_harm
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_tor_harm SHARED dtorh1.f dtorh2.f dtorh3.f rout.f)

add_executable(test_tor_harm testpro.f)
target_link_libraries(test_tor_harm gst_tor_harm)
