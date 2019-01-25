cmake_minimum_required (VERSION 3.0.2)

project(
  JacobiTest
  VERSION 1.0.0
  LANGUAGES C CXX
)

set(JACOBI_ROURCES src/jacobi.cpp src/gauss_quad.cpp src/interface.cpp)

add_executable(deriv test/deriv.cpp ${JACOBI_ROURCES})
target_link_libraries(deriv gsl blas)

add_executable(integrate test/integrate.cpp ${JACOBI_ROURCES})
target_link_libraries(integrate gsl blas)

add_executable(interp test/interp.cpp ${JACOBI_ROURCES})
target_link_libraries(interp gsl blas)
