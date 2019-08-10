cmake_minimum_required (VERSION 3.10)

project(
  toms549_weierstrass_p
  VERSION 1.0
  LANGUAGES Fortran
)

enable_testing()

add_library(toms549_weierstrass_p SHARED toms549.f)

add_executable(test_toms549_equiharmonic test_equiharmonic.f)
target_link_libraries(test_toms549_equiharmonic toms549_weierstrass_p)

add_executable(test_toms549_lemniscatic test_lemniscatic.f)
target_link_libraries(test_toms549_lemniscatic toms549_weierstrass_p)
