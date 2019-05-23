cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

add_executable(test_appell_f1 test_appell_f1.cpp)
#target_include_directories(test_appell_f1 PRIVATE )
#target_link_libraries( )
