cmake_minimum_required (VERSION 3.10)

project(
  gst_pro_obl_harm
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(gst_pro_obl_harm SHARED doblh.f dproh.f fact.f sph_harm.f95)

add_executable(test_pro_obl_harm testpro.f)
target_link_libraries(test_pro_obl_harm gst_pro_obl_harm)

add_executable(test_sph_harm test_sph_harm.cpp)
target_link_libraries(test_sph_harm gst_pro_obl_harm)
