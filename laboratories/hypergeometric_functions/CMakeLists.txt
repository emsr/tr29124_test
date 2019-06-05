cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_hypergeometric_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_CXX_FLAGS -Wno-psabi)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing( )

add_executable(hyperg_1F2 hyperg_1F2.cpp)
target_compile_definitions(hyperg_1F2 PRIVATE MAIN_1F2)
#target_include_directories(hyperg_1F2 PRIVATE )
#target_link_libraries(hyperg_1F2 )

add_executable(hyperg_2F0 hyperg_2F0.cpp)
target_compile_definitions(hyperg_2F0 PRIVATE MAIN_2F0)
#target_include_directories(hyperg_2F0 PRIVATE )
#target_link_libraries(hyperg_2F0 )

add_executable(hyperg_3F0 hyperg_3F0.cpp)
target_compile_definitions(hyperg_3F0 PRIVATE MAIN_3F0)
#target_include_directories(hyperg_3F0 PRIVATE )
#target_link_libraries(hyperg_3F0 )

add_executable(test_conf_hyperg test_conf_hyperg.cpp)
target_include_directories(test_conf_hyperg PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_conf_hyperg quadmath)

add_executable(test_conf_hyperg_limit test_conf_hyperg_limit.cpp)
target_include_directories(test_conf_hyperg_limit PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_conf_hyperg_limit quadmath)

add_executable(test_hyperg test_hyperg.cpp)
target_include_directories(test_hyperg PRIVATE
  ../../wrappers
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_hyperg
  quadmath
  -L../../wrappers/debug -lwrap_gsl)

#add_executable(test_tricomi_u test_tricomi_u.cpp)
#target_include_directories(test_tricomi_u PRIVATE )
#target_link_libraries(test_tricomi_u )
