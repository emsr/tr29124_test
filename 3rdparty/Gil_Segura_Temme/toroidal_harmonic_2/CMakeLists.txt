cmake_minimum_required (VERSION 3.10)

project(
  gst_toroidal_harmonic_2
  VERSION 2.0
  LANGUAGES Fortran
)

add_library(gst_toroidal_harmonic_2 SHARED adrout.f dtorh3_2.f Rout.f)
#target_include_directories(gst_toroidal_harmonic_2 PUBLIC .)

