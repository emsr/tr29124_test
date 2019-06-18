cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_airy_functions
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

enable_testing()

add_executable(airy_toy airy_toy.cpp)
target_include_directories(airy_toy PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(airy_toy quadmath)
add_test(NAME run_airy_toy COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/airy_toy > output/airy_toy.txt")

add_executable(airy_toy_old airy_toy_old.cpp)
target_include_directories(airy_toy_old PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(airy_toy_old quadmath)
add_test(NAME run_airy_toy_old COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/airy_toy_old > output/airy_toy_old.txt")

add_executable(debug_airy debug_airy.cpp)
target_include_directories(debug_airy PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(debug_airy quadmath -L../../wrappers -lwrap_gsl)
add_test(NAME run_debug_airy COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/debug_airy > output/debug_airy.txt")

add_executable(plot_airy plot_airy.cpp)
target_include_directories(plot_airy PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(plot_airy quadmath)
add_test(NAME run_plot_airy COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/plot_airy > output/plot_airy.txt")

add_executable(test_airy_roots test_airy_roots.cpp)
target_include_directories(test_airy_roots PRIVATE
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
)
target_link_libraries(test_airy_roots quadmath)
add_test(NAME run_test_airy_roots COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_airy_roots > output/test_airy_roots.txt")

add_executable(test_little_airy test_little_airy.cpp)
target_include_directories(test_little_airy PRIVATE
  .
  ../../cxx_math_constants/include
  ../../include
  ../../cxx_fp_utils/include
  ../../cxx_complex_utils/include
  ../../cxx_complex_math/include
  ../../cxx_summation/include
  ../../polynomial/include
  ../../quadrature/include
  ../../wrappers
)
target_link_libraries(test_little_airy quadmath -L../../wrappers -lwrap_gsl)
add_test(NAME run_test_little_airy COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_little_airy > output/test_little_airy.txt")
