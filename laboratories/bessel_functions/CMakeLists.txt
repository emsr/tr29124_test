cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_bessel_functions
  VERSION 0.5
  LANGUAGES CXX Fortran
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_custom_target(make_lab_bessel_functions_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(debug_bessel debug_bessel.cpp)

add_executable(test_cyl_bessel_ratio test_cyl_bessel_ratio.cpp)
target_link_libraries(test_cyl_bessel_ratio ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_cyl_bessel_ratio COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cyl_bessel_ratio > output/test_cyl_bessel_ratio.txt")

add_executable(test_cyl_hankel_ratio test_cyl_hankel_ratio.cpp)
target_link_libraries(test_cyl_hankel_ratio ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_cyl_hankel_ratio COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cyl_hankel_ratio > output/test_cyl_hankel_ratio.txt")

add_executable(hankel_toy128 hankel_toy128.cpp)
target_link_libraries(hankel_toy128 ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_hankel_toy128 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hankel_toy128 > output/hankel_toy128.txt")

add_executable(hankel_toy hankel_toy.cpp)
target_link_libraries(hankel_toy ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_hankel_toy COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hankel_toy > output/hankel_toy.txt")

add_executable(hankel_toy_new hankel_toy_new.cpp)
target_link_libraries(hankel_toy_new ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_hankel_toy_new COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hankel_toy_new > output/hankel_toy_new.txt")

#add_executable(hankel_transition hankel_transition.cpp)
#target_link_libraries(hankel_transition ${EMSR_LIBRARIES} quadmath)
#add_test(NAME run_hankel_transition COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/hankel_transition > output/hankel_transition.txt")

add_executable(help_kelvin help_kelvin.cpp)
target_link_libraries(help_kelvin ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_help_kelvin COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/help_kelvin > output/help_kelvin.txt")

add_executable(plot_bessel plot_bessel.cpp)
target_link_libraries(plot_bessel ${EMSR_LIBRARIES})
add_test(NAME run_plot_bessel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/plot_bessel > output/plot_bessel.txt")

add_executable(test_anger_weber test_anger_weber.cpp)
target_link_libraries(test_anger_weber ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_anger_weber COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_anger_weber > output/test_anger_weber.txt")

add_executable(test_bessel_asymp test_bessel_asymp.cpp)
target_link_libraries(test_bessel_asymp ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_bessel_asymp COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bessel_asymp > output/test_bessel_asymp.txt")

add_executable(test_bessel test_bessel.cpp)
target_link_libraries(test_bessel ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_bessel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bessel > output/test_bessel.txt")

add_executable(test_bessel_iter test_bessel_iter.cpp)
target_link_libraries(test_bessel_iter ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_bessel_iter COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bessel_iter > output/test_bessel_iter.txt")

add_executable(test_cyl_bessel_zeros test_cyl_bessel_zeros.cpp)
target_link_libraries(test_cyl_bessel_zeros ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_cyl_bessel_zeros COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cyl_bessel_zeros > output/test_cyl_bessel_zeros.txt")

add_executable(test_hankel test_hankel.cpp)
target_link_libraries(test_hankel ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_hankel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hankel > output/test_hankel.txt")

if (Boost_FOUND)
  add_executable(test_cyl_hankel test_cyl_hankel.cpp)
  target_link_libraries(test_cyl_hankel quadmath ${EMSR_LIBRARIES} wrap_boost)
  add_test(NAME run_test_cyl_hankel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cyl_hankel > output/test_cyl_hankel.txt")

  add_executable(test_hankel_real_arg test_hankel_real_arg.cpp)
  target_link_libraries(test_hankel_real_arg ${EMSR_LIBRARIES} quadmath wrap_boost)
  add_test(NAME run_test_hankel_real_arg COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hankel_real_arg > output/test_hankel_real_arg.txt")
endif(Boost_FOUND)

add_executable(test_kelvin test_kelvin.cpp)
target_link_libraries(test_kelvin ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_kelvin COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_kelvin > output/test_kelvin.txt")

add_executable(test_large_order_bessel test_large_order_bessel.cpp)
target_link_libraries(test_large_order_bessel ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_large_order_bessel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_large_order_bessel > output/test_large_order_bessel.txt")

add_executable(test_lommel test_lommel.cpp)
target_link_libraries(test_lommel ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_lommel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lommel > output/test_lommel.txt")

add_executable(test_mod_bessel_asymp test_mod_bessel_asymp.cpp)
target_link_libraries(test_mod_bessel_asymp ${EMSR_LIBRARIES} cxx_float128 quadmath)

add_executable(test_nric_bessel test_nric_bessel.cpp)
target_include_directories(test_nric_bessel PRIVATE
  .)
add_test(NAME run_test_nric_bessel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_nric_bessel > output/test_nric_bessel.txt")

add_executable(test_sph_bessel test_sph_bessel.cpp)
target_link_libraries(test_sph_bessel ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_sph_bessel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sph_bessel > output/test_sph_bessel.txt")

if (Boost_FOUND)
  add_executable(test_sph_hankel test_sph_hankel.cpp)
  target_link_libraries(test_sph_hankel ${EMSR_LIBRARIES} cxx_float128
    quadmath
    wrap_boost)
  add_test(NAME run_test_sph_hankel COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_sph_hankel > output/test_sph_hankel.txt")
endif(Boost_FOUND)

add_executable(test_struve test_struve.cpp)
target_include_directories(test_struve PRIVATE
  .)
target_link_libraries(test_struve quadmath ${EMSR_LIBRARIES} cxx_float128 wrap_burkhardt gfortran)
#FIXME add_test(NAME run_test_struve COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_struve > output/test_struve.txt")

add_executable(test_cyl_bessel_nk_series test_cyl_bessel_nk_series.cpp)
target_link_libraries(test_cyl_bessel_nk_series ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_cyl_bessel_nk_series COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_cyl_bessel_nk_series > output/test_cyl_bessel_nk_series.txt")

add_executable(test_struve_old test_struve_old.cpp)
target_include_directories(test_struve_old PRIVATE
  .
  ../hypergeometric_functions)
target_link_libraries(test_struve_old ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_test_struve_old COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_struve_old > output/test_struve_old.txt")
