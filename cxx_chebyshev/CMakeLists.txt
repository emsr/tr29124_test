cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_chebyshev
  VERSION 0.1
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

add_custom_target(make_cxx_chebyshev_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_library(cxx_chebyshev INTERFACE)
target_include_directories(cxx_chebyshev INTERFACE include)

add_executable(test_cxx_chebyshev test_chebyshev.cpp)
target_link_libraries(test_cxx_chebyshev
  cxx_math_constants
  cxx_float128
  cxx_complex_utils
  cxx_polynomial
  cxx_fp_utils
  cxx_chebyshev
  quadmath)

add_test(NAME run_test_cxx_chebyshev COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cxx_chebyshev > output/test_cxx_chebyshev.txt")

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_cxx_chebyshev
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_chebyshev" VERBATIM
  )
endif (DOXYGEN_FOUND)
