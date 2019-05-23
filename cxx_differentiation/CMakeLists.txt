cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_differentiation
  VERSION 1.0.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

list(APPEND cxx_differentiation_INCLUDES
  include/ext/derivative.h
  include/ext/derivative.tcc
)


add_executable(test_derivative test_derivative.cpp)
target_include_directories(test_derivative PRIVATE include)

add_test(NAME derivative_test0 COMMAND bash -c "./test_derivative > test_derivative.txt")

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_differentiation" VERBATIM
  )
endif (DOXYGEN_FOUND)
