cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_elliptic_integrals
  VERSION 0.5
  LANGUAGES CXX
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

add_custom_target(make_lab_elliptic_integrals_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

#add_executable(carlson1 carlson1.cpp)
#target_include_directories(carlson1 PRIVATE )
#target_link_libraries(carlson1 )

#add_executable(carlson2 carlson2.cpp)
#target_include_directories(carlson2 PRIVATE )
#target_link_libraries(carlson2 )

#add_executable(carlson3 carlson3.cpp)
#target_include_directories(carlson3 PRIVATE )
#target_link_libraries(carlson3 )

#add_executable(carlson4 carlson4.cpp)
#target_include_directories(carlson4 PRIVATE )
#target_link_libraries(carlson4 )

#add_executable(carlson carlson.cpp)
#target_include_directories(carlson PRIVATE )
#target_link_libraries(carlson )

add_executable(carlson_std carlson_std.cpp)
target_link_libraries(carlson_std ${EMSR_LIBRARIES} quadmath)
add_test(NAME run_carlson_std COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/carlson_std > output/carlson_std.txt")

add_executable(test_comp_ellint_1 test_comp_ellint_1.cpp)
target_link_libraries(test_comp_ellint_1 ${EMSR_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_comp_ellint_1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_comp_ellint_1 > output/test_comp_ellint_1.txt")

if (Boost_FOUND)
  add_executable(test_heuman_lambda test_heuman_lambda.cpp)
  target_link_libraries(test_heuman_lambda ${EMSR_LIBRARIES} quadmath wrap_boost)
  add_test(NAME run_test_heuman_lambda COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_heuman_lambda > output/test_heuman_lambda.txt")

  add_executable(test_jacobi_zeta test_jacobi_zeta.cpp)
  target_link_libraries(test_jacobi_zeta ${EMSR_LIBRARIES} quadmath wrap_boost)
  add_test(NAME run_test_jacobi_zeta COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_zeta > output/test_jacobi_zeta.txt")
endif(Boost_FOUND)

add_executable(test_legendre_ellint test_legendre_ellint.cpp)
target_link_libraries(test_legendre_ellint ${EMSR_LIBRARIES} quadmath wrap_gsl)
add_test(NAME run_test_legendre_ellint COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre_ellint > output/test_legendre_ellint.txt")
