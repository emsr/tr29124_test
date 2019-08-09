cmake_minimum_required (VERSION 3.10)

project(
  gst_tor_harm2
  VERSION 2.0
  LANGUAGES Fortran
)

add_library(gst_tor_harm2 SHARED adrout.f dtorh3_2.f Rout.f)

