cmake_minimum_required (VERSION 3.0.2)

project(
  QuarticTest
  VERSION 1.0.0
  LANGUAGES C CXX
)

add_executable(quartic test.c solve_quartic.c zsolve_quartic.c)
target_link_libraries(quartic gsl blas m)
