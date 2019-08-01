cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_continued_fractions
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_library(cxx_continued_fractions INTERFACE)
target_include_directories(cxx_continued_fractions INTERFACE include)

add_executable(test_continued_fraction test_continued_fraction.cpp)
target_include_directories(test_continued_fraction PRIVATE
  ../include # For bits/numeric_limits.h
  ../cxx_float128/include)
target_link_libraries(test_continued_fraction cxx_continued_fractions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_summation
  cxx_polynomial
  cxx_integration
  cxx_special_functions
  quadmath)
add_test(NAME run_test_continued_fraction COMMAND bash -c "./bin/test_continued_fraction > test_continued_fraction.txt")

add_executable(test_lentz_continued_fraction test_lentz_continued_fraction.cpp)
target_include_directories(test_lentz_continued_fraction PRIVATE
  ../include # For bits/numeric_limits.h
)
target_link_libraries(test_lentz_continued_fraction
  cxx_continued_fractions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_summation
  cxx_polynomial
  cxx_integration
  cxx_special_functions
  quadmath)
add_test(NAME run_test_lentz_continued_fraction COMMAND bash -c "./bin/test_lentz_continued_fraction > test_lentz_continued_fraction.txt")

add_executable(test_steed_continued_fraction test_steed_continued_fraction.cpp)
target_include_directories(test_steed_continued_fraction PRIVATE
  ../include # For bits/numeric_limits.h
)
target_link_libraries(test_steed_continued_fraction
  cxx_continued_fractions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_summation
  cxx_polynomial
  cxx_integration
  cxx_special_functions
  quadmath)
add_test(NAME run_test_steed_continued_fraction COMMAND bash -c "./bin/test_steed_continued_fraction > test_steed_continued_fraction.txt")

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_continued_fractions" VERBATIM
  )
endif (DOXYGEN_FOUND)
