cmake_minimum_required (VERSION 3.10)

project(
  gst_toroidal_harmonic
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_toroidal_harmonic SHARED dtorh1.f dtorh2.f dtorh3.f rout.f toroidal_harmonic.f90)
target_include_directories(gst_toroidal_harmonic PUBLIC .)

add_executable(test_tor_harm testpro.f)
target_link_libraries(test_tor_harm gst_toroidal_harmonic)

add_executable(test_gst_toroidal_harmonic test_toroidal_harmonic.cpp)
target_link_libraries(test_gst_toroidal_harmonic gst_toroidal_harmonic)
