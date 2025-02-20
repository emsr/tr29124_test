cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_exponential_integrals
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_lab_exponential_integrals_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(test_csint test_csint.cpp)
target_link_libraries(test_csint ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_csint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_csint > output/test_csint.txt")

if (HAVE_BOOST)
  add_executable(test_expint test_expint.cpp)
  target_link_libraries(test_expint
    ${EMSR_LIBRARIES}
    cxx_float128
    quadmath
    wrap_boost)
  add_test(NAME run_test_expint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_expint > output/test_expint.txt")
endif(HAVE_BOOST)
