cmake_minimum_required (VERSION 3.10)

project(
  toms804_adapt_dubdiv_quad
  VERSION 1.0
  LANGUAGES Fortran
)

enable_testing()

set(lib_src_dir F90/Src)
set(test_src_dir F90/Drivers)

add_library(toms954_low_degree_solvers SHARED ${lib_src_dir}/Polynomial234RootSolvers.f90 ${lib_src_dir}/SetWorkingPrecision.f90)

add_executable(comp_toms954_low_degree_solvers ${test_src_dir}/PolynomialRootsCompare.f90 ${test_src_dir}/JenkinsTraubSolver.f90)
target_link_libraries(comp_toms954_low_degree_solvers toms954_low_degree_solvers)

add_executable(time_toms954_low_degree_solvers ${test_src_dir}/PolynomialRootsTimings.f90 ${test_src_dir}/JenkinsTraubSolver.f90)
target_link_libraries(time_toms954_low_degree_solvers toms954_low_degree_solvers)
