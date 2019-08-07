cmake_minimum_required (VERSION 3.10)

project(
  incgam
  VERSION 1.0
  LANGUAGES Fortran
)

add_library(incgam SHARED cons.f90   incgam.f90)

add_executable(testincgam testincgam.f90)
target_link_libraries(testincgam incgam)
