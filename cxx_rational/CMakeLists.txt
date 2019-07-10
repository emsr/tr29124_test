cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_rational
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_library(cxx_rational INTERFACE)
target_include_directories(cxx_rational INTERFACE include)

add_executable(test_rational test_rational.cpp)
target_link_libraries(test_rational cxx_rational quadmath)
target_include_directories(test_rational PRIVATE
  ../include # float128
  ../cxx_complex_utils/include
  ../cxx_polynomial/include
  ../cxx_fp_utils/include)

add_test(NAME run_test_rational COMMAND bash -c "./bin/test_rational > output/test_rational.txt")

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_rational
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_rational" VERBATIM
  )
endif (DOXYGEN_FOUND)
