cmake_minimum_required (VERSION 3.10)

project(
  gst_scorer
  VERSION 1.0
  LANGUAGES Fortran
)

set(src_dir Fortran77/Src/Dp)
add_library(gst_scorer SHARED ${src_dir}/src_819.f ${src_dir}/src.f)
