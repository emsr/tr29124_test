cmake_minimum_required (VERSION 3.10.0)

project(
  cxx_chebyshev
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

#add_library( INTERFACE)
#target_include_directories( INTERFACE qwerty)

add_executable(testcase testcase.cpp)
target_link_libraries(testcase -L../wrappers quadmath -lwrap_gsl -lwrap_boost -lwrap_burkhardt)
target_include_directories(testcase PRIVATE
  .
  ../wrappers
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include
  ../quadrature/include
  ../cxx_math_const/include)

add_executable(testcase_tr1 testcase.cpp)
target_link_libraries(testcase_tr1 quadmath -L../wrappers -lwrap_gsl -lwrap_boost)
target_compile_definitions(testcase_tr1 PRIVATE TR1)
target_include_directories(testcase_tr1 PRIVATE
  .
  wrappers
  ../include
  ../cxx_fp_utils/include
  ../polynomial/include
  ../cxx_summation/include
  ../quadrature/include
  ../cxx_math_const/include)
