cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_bessel_functions
  VERSION 0.5
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "/home/ed/bin/bin/g++")
set(CMAKE_FORTRAN_COMPILER "/home/ed/bin/bin/gfortran")
set(CMAKE_CXX_FLAGS -Wno-psabi)

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

enable_testing()

add_executable(debug_bessel debug_bessel.cpp)

add_executable(test_cyl_bessel_ratio test_cyl_bessel_ratio.cpp)
target_include_directories(test_cyl_bessel_ratio PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_continued_fractions/include)

add_executable(test_cyl_hankel_ratio test_cyl_hankel_ratio.cpp)
target_include_directories(test_cyl_hankel_ratio PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_continued_fractions/include)

add_executable(hankel_toy128 hankel_toy128.cpp)
target_include_directories(hankel_toy128 PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include)
target_link_libraries(hankel_toy128 quadmath)

add_executable(hankel_toy hankel_toy.cpp)
target_include_directories(hankel_toy PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include)
target_link_libraries(hankel_toy quadmath)

add_executable(hankel_toy_new hankel_toy_new.cpp)
target_include_directories(hankel_toy_new PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include)
target_link_libraries(hankel_toy_new quadmath)

#add_executable(hankel_transition hankel_transition.cpp)
#target_include_directories(hankel_transition PRIVATE
#  ../../include
#  ../../cxx_fp_utils/include
#  ../../cxx_summation/include
#  ../../polynomial/include
#  ../../quadrature/include)
#target_link_libraries(hankel_transition quadmath)

add_executable(help_kelvin help_kelvin.cpp)
target_include_directories(help_kelvin PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(help_kelvin quadmath)

add_executable(plot_bessel plot_bessel.cpp)

add_executable(test_anger_weber test_anger_weber.cpp)
target_include_directories(test_anger_weber PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_anger_weber quadmath)

add_executable(test_bessel_asymp test_bessel_asymp.cpp)
target_include_directories(test_bessel_asymp PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_bessel_asymp quadmath)

add_executable(test_bessel test_bessel.cpp)
target_include_directories(test_bessel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_bessel quadmath)

add_executable(test_bessel_iter test_bessel_iter.cpp)
target_include_directories(test_bessel_iter PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_bessel_iter quadmath)

add_executable(test_cyl_bessel_zeros test_cyl_bessel_zeros.cpp)
target_include_directories(test_cyl_bessel_zeros PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_cyl_bessel_zeros quadmath)

add_executable(test_cyl_hankel test_cyl_hankel.cpp)
target_include_directories(test_cyl_hankel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_cyl_hankel quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_hankel test_hankel.cpp)
target_include_directories(test_hankel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_hankel quadmath)

add_executable(test_hankel_real_arg test_hankel_real_arg.cpp)
target_include_directories(test_hankel_real_arg PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_hankel_real_arg quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_kelvin test_kelvin.cpp)
target_include_directories(test_kelvin PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_kelvin quadmath)

add_executable(test_large_order_bessel test_large_order_bessel.cpp)
target_include_directories(test_large_order_bessel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_large_order_bessel quadmath)

add_executable(test_lommel test_lommel.cpp)
target_include_directories(test_lommel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_lommel quadmath)

add_executable(test_mod_bessel_asymp test_mod_bessel_asymp.cpp)
target_include_directories(test_mod_bessel_asymp PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_mod_bessel_asymp quadmath)

add_executable(test_nric_bessel test_nric_bessel.cpp)
target_include_directories(test_nric_bessel PRIVATE
  .
)

add_executable(test_sph_bessel test_sph_bessel.cpp)
target_include_directories(test_sph_bessel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_sph_bessel quadmath)

add_executable(test_sph_hankel test_sph_hankel.cpp)
target_include_directories(test_sph_hankel PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include
  ../../wrappers)
target_link_libraries(test_sph_hankel quadmath -L../../wrappers/debug -lwrap_boost)

add_executable(test_struve test_struve.cpp)
target_include_directories(test_struve PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include
  ../../wrappers)
target_link_libraries(test_struve quadmath -L../../wrappers/debug -lwrap_burkhardt -L/home/ed/bin/lib64 -lgfortran)

add_executable(test_cyl_bessel_nk_series test_cyl_bessel_nk_series.cpp)
target_include_directories(test_cyl_bessel_nk_series PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../quadrature/include
  ../../cxx_summation/include)
target_link_libraries(test_cyl_bessel_nk_series quadmath)

# Missing hyperg_1F2.cpp
#add_executable(test_struve_old test_struve_old.cpp)
#target_include_directories(test_struve_old PRIVATE )
#target_link_libraries(test_struve_old )
