cmake_minimum_required (VERSION 3.10)

project(
  gst_spheroidal_harmonic
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_spheroidal_harmonic SHARED doblh.f dproh.f fact.f sph_harm.f95)
target_include_directories(gst_spheroidal_harmonic PUBLIC .)

add_executable(test_gst_spheroidal_harmonic testpro.f)
target_link_libraries(test_gst_spheroidal_harmonic gst_spheroidal_harmonic)

add_executable(test_sph_harm test_sph_harm.cpp)
target_link_libraries(test_sph_harm gst_spheroidal_harmonic)
