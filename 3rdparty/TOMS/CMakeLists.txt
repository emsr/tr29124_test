cmake_minimum_required (VERSION 3.10)

project(
  toms
  VERSION 1.0
  LANGUAGES CXX
)

add_subdirectory(toms353_filon_quad)
add_subdirectory(toms468_adapt_subdiv_quad)
add_subdirectory(toms549_weierstrass_p)
add_subdirectory(toms804_mathieu)
add_subdirectory(toms812_bernstein)
add_subdirectory(toms917_wright_omega)
