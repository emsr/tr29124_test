cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_differentiation
  VERSION 1.0.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_library(cxx_differentiation INTERFACE)
target_include_directories(cxx_differentiation INTERFACE include)

add_executable(test_differentiation test_differentiation.cpp)
target_link_libraries(test_differentiation cxx_differentiation)

add_test(NAME run_test_differentiation COMMAND bash -c "./bin/test_differentiation > output/test_differentiation.txt")

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_differentiation" VERBATIM
  )
endif (DOXYGEN_FOUND)
