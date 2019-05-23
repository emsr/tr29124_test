cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  specfun_laboratories
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

add_executable(test_q_number test_q_number.cpp)
#target_include_directories(test_q_number PRIVATE )
#target_link_libraries(test_q_number )

add_executable(test_q_pochhammer test_q_pochhammer.cpp)
#target_include_directories(test_q_pochhammer PRIVATE )
#target_link_libraries(test_q_pochhammer )
