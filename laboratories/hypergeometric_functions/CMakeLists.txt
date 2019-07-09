cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_hypergeometric_functions
  VERSION 0.5
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

enable_testing( )

add_executable(hyperg_1F2 hyperg_1F2.cpp)
target_compile_definitions(hyperg_1F2 PRIVATE MAIN_1F2)
add_test(NAME run_hyperg_1F2 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hyperg_1F2 > output/hyperg_1F2.txt")

add_executable(hyperg_2F0 hyperg_2F0.cpp)
target_compile_definitions(hyperg_2F0 PRIVATE MAIN_2F0)
add_test(NAME run_hyperg_2F0 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hyperg_2F0 > output/hyperg_2F0.txt")

add_executable(hyperg_3F0 hyperg_3F0.cpp)
target_compile_definitions(hyperg_3F0 PRIVATE MAIN_3F0)
add_test(NAME run_hyperg_3F0 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hyperg_3F0 > output/hyperg_3F0.txt")

add_executable(test_conf_hyperg test_conf_hyperg.cpp)
target_include_directories(test_conf_hyperg PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_conf_hyperg quadmath)
add_test(NAME run_test_conf_hyperg COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_conf_hyperg > output/test_conf_hyperg.txt")

add_executable(test_conf_hyperg_limit test_conf_hyperg_limit.cpp)
target_include_directories(test_conf_hyperg_limit PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_conf_hyperg_limit quadmath)
add_test(NAME run_test_conf_hyperg_limit COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_conf_hyperg_limit > output/test_conf_hyperg_limit.txt")

add_executable(test_hyperg test_hyperg.cpp)
target_include_directories(test_hyperg PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include
  ../../wrappers)
target_link_libraries(test_hyperg
  quadmath
  -L../../wrappers -lwrap_gsl)
add_test(NAME run_test_hyperg COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hyperg > output/test_hyperg.txt")

#add_executable(test_tricomi_u test_tricomi_u.cpp)
#target_include_directories(test_tricomi_u PRIVATE )
#target_link_libraries(test_tricomi_u )
