cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_continued_fractions
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_executable(test_continued_fraction test_continued_fraction.cpp)
target_include_directories(test_continued_fraction PRIVATE
  include
  ../include
  ../cxx_fp_utils/include
  ../cxx_summation/include
  ../polynomial/include
  ../quadrature/include)
target_link_libraries(test_continued_fraction quadmath)

add_executable(test_lentz_continued_fraction test_lentz_continued_fraction.cpp)
target_include_directories(test_lentz_continued_fraction PRIVATE
  include
  ../include
  ../cxx_fp_utils/include
  ../cxx_summation/include
  ../polynomial/include
  ../quadrature/include)
target_link_libraries(test_lentz_continued_fraction quadmath)

add_executable(test_steed_continued_fraction test_steed_continued_fraction.cpp)
target_include_directories(test_steed_continued_fraction PRIVATE
  include
  ../include
  ../cxx_fp_utils/include
  ../cxx_summation/include
  ../polynomial/include
  ../quadrature/include)
target_link_libraries(test_steed_continued_fraction quadmath)
