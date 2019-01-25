cmake_minimum_required (VERSION 3.0.2)

project(
  FresnelTest
  VERSION 1.0.0
  LANGUAGES C CXX
)

add_executable(test_fresnel test_fresnel.c fresnel.c)
target_link_libraries(test_fresnel gsl blas m)
