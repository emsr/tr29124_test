cmake_minimum_required (VERSION 3.10)

project(
  toms812_bernstein
  VERSION 1.0
  LANGUAGES CXX
)

enable_testing()

set(lib_src_dir Cpp/Sp/Src)
set(test_src_dir Cpp/Sp/Drivers)

add_library(toms812_bernstein SHARED ${lib_src_dir}/Bernstein.cpp)
target_include_directories(toms812_bernstein PUBLIC ${lib_src_dir})

add_executable(test_toms812_bernstein_1 ${test_src_dir}/driver1.cpp)
target_link_libraries(test_toms812_bernstein_1 toms812_bernstein)

add_executable(test_toms812_bernstein_2 ${test_src_dir}/driver2.cpp)
target_link_libraries(test_toms812_bernstein_2 toms812_bernstein)

add_executable(test_toms812_bernstein_3 ${test_src_dir}/driver3.cpp)
target_link_libraries(test_toms812_bernstein_3 toms812_bernstein)

add_executable(test_toms812_bernstein_4 ${test_src_dir}/driver4.cpp)
target_link_libraries(test_toms812_bernstein_4 toms812_bernstein)

add_executable(test_toms812_bernstein_5 ${test_src_dir}/driver5.cpp)
target_link_libraries(test_toms812_bernstein_5 toms812_bernstein)

add_executable(test_toms812_bernstein_6 ${test_src_dir}/driver6.cpp)
target_link_libraries(test_toms812_bernstein_6 toms812_bernstein)

add_executable(test_toms812_bernstein_7 ${test_src_dir}/driver7.cpp)
target_link_libraries(test_toms812_bernstein_7 toms812_bernstein)

add_executable(test_toms812_bernstein_8 ${test_src_dir}/driver8.cpp)
target_link_libraries(test_toms812_bernstein_8 toms812_bernstein)

add_executable(test_toms812_bernstein_9 ${test_src_dir}/driver9.cpp)
target_link_libraries(test_toms812_bernstein_9 toms812_bernstein)

add_executable(test_toms812_bernstein_10 ${test_src_dir}/driver10.cpp)
target_link_libraries(test_toms812_bernstein_10 toms812_bernstein)
