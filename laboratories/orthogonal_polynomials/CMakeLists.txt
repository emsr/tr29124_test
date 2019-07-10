cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_orthogonal_polynomials
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

add_executable(debug_odd_legendre debug_odd_legendre.cpp)
target_include_directories(debug_odd_legendre PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(debug_odd_legendre quadmath)
add_test(NAME run_debug_odd_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/debug_odd_legendre > output/debug_odd_legendre.txt")

add_executable(test_assoc_laguerre test_assoc_laguerre.cpp)
target_include_directories(test_assoc_laguerre PRIVATE )
target_link_libraries(test_assoc_laguerre quadmath)
add_test(NAME run_test_assoc_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_laguerre > output/test_assoc_laguerre.txt")

add_executable(test_assoc_legendre test_assoc_legendre.cpp)
target_include_directories(test_assoc_legendre PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_assoc_legendre quadmath)
add_test(NAME run_test_assoc_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_legendre > output/test_assoc_legendre.txt")

add_executable(test_assoc_legendre_q test_assoc_legendre_q.cpp)
target_include_directories(test_assoc_legendre_q PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_assoc_legendre_q quadmath)
add_test(NAME run_test_assoc_legendre_q COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_legendre_q > output/test_assoc_legendre_q.txt")

add_executable(test_charlier test_charlier.cpp)
target_include_directories(test_charlier PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../wrappers)
target_link_libraries(test_charlier gfortran quadmath -L../../wrappers -lwrap_burkhardt)
add_test(NAME run_test_charlier COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_charlier >test_charlier .txt")

add_executable(test_chebyshev test_chebyshev.cpp)
target_include_directories(test_chebyshev PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_chebyshev quadmath)
add_test(NAME run_test_chebyshev COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev > output/test_chebyshev.txt")

add_executable(test_chebyshev_trig test_chebyshev_trig.cpp)
target_include_directories(test_chebyshev_trig PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_chebyshev_trig quadmath)
add_test(NAME run_test_chebyshev_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev_trig > output/test_chebyshev_trig.txt")

add_executable(test_chebyshev_trig_pi test_chebyshev_trig_pi.cpp)
target_include_directories(test_chebyshev_trig_pi PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_chebyshev_trig_pi quadmath)
add_test(NAME run_test_chebyshev_trig_pi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev_trig_pi > output/test_chebyshev_trig_pi.txt")

add_executable(test_complex test_complex.cpp)
target_include_directories(test_complex PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_complex quadmath)
add_test(NAME run_test_complex COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_complex > output/test_complex.txt")

add_executable(test_continuous_dual_hahn test_continuous_dual_hahn.cpp)
target_include_directories(test_continuous_dual_hahn PRIVATE )
target_link_libraries(test_continuous_dual_hahn quadmath)
add_test(NAME run_test_continuous_dual_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_continuous_dual_hahn > output/test_continuous_dual_hahn.txt")

add_executable(test_continuous_hahn test_continuous_hahn.cpp)
target_include_directories(test_continuous_hahn PRIVATE )
target_link_libraries(test_continuous_hahn quadmath)

add_executable(test_dual_hahn test_dual_hahn.cpp)
target_include_directories(test_dual_hahn PRIVATE )
target_link_libraries(test_dual_hahn quadmath)
add_test(NAME run_test_dual_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dual_hahn > output/test_dual_hahn.txt")

add_executable(test_gegenbauer test_gegenbauer.cpp)
target_include_directories(test_gegenbauer PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_gegenbauer quadmath)
add_test(NAME run_test_gegenbauer COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer > output/test_gegenbauer.txt")

add_executable(test_gegenbauer_neg_parm test_gegenbauer_neg_parm.cpp)
target_include_directories(test_gegenbauer_neg_parm PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_gegenbauer_neg_parm quadmath)
add_test(NAME run_test_gegenbauer_neg_parm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer_neg_parm > output/test_gegenbauer_neg_parm.txt")

add_executable(test_gegenbauer_neg_roots test_gegenbauer_neg_roots.cpp)
target_include_directories(test_gegenbauer_neg_roots PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_gegenbauer_neg_roots quadmath)
add_test(NAME run_test_gegenbauer_neg_roots COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer_neg_roots > output/test_gegenbauer_neg_roots.txt")

add_executable(test_hahn test_hahn.cpp)
target_include_directories(test_hahn PRIVATE )
target_link_libraries(test_hahn quadmath)
add_test(NAME run_test_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hahn > output/test_hahn.txt")

add_executable(test_hermite test_hermite.cpp)
target_include_directories(test_hermite PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_hermite quadmath)
add_test(NAME run_test_hermite COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hermite > output/test_hermite.txt")

add_executable(test_jacobi test_jacobi.cpp)
target_include_directories(test_jacobi PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_jacobi quadmath)
add_test(NAME run_test_jacobi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi > output/test_jacobi.txt")

add_executable(test_jacobi_neg_parm test_jacobi_neg_parm.cpp)
target_include_directories(test_jacobi_neg_parm PRIVATE
  .
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_jacobi_neg_parm quadmath)
add_test(NAME run_test_jacobi_neg_parm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_neg_parm > output/test_jacobi_neg_parm.txt")

add_executable(test_jacobi_neg_roots test_jacobi_neg_roots.cpp)
target_include_directories(test_jacobi_neg_roots PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_jacobi_neg_roots quadmath)
add_test(NAME run_test_jacobi_neg_roots COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_neg_roots > output/test_jacobi_neg_roots.txt")

add_executable(test_krawtchouk test_krawtchouk.cpp)
target_include_directories(test_krawtchouk PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../wrappers)
target_link_libraries(test_krawtchouk
  gfortran
  quadmath
  -L../../wrappers -lwrap_burkhardt)
add_test(NAME run_test_krawtchouk COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_krawtchouk > output/test_krawtchouk.txt")

add_executable(test_laguerre test_laguerre.cpp)
target_include_directories(test_laguerre PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_laguerre quadmath)
add_test(NAME run_test_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_laguerre > output/test_laguerre.txt")

add_executable(test_legendre test_legendre.cpp)
target_include_directories(test_legendre PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_legendre quadmath)
add_test(NAME run_test_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre > output/test_legendre.txt")

add_executable(test_legendre_q test_legendre_q.cpp)
target_include_directories(test_legendre_q PRIVATE
  ${TR29124_INCLUDE_DIRS})
target_link_libraries(test_legendre_q quadmath)
add_test(NAME run_test_legendre_q COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre_q > output/test_legendre_q.txt")

add_executable(test_lobatto test_lobatto.cpp)
target_include_directories(test_lobatto PRIVATE
  .
  ${TR29124_INCLUDE_DIRS}
  ../../cxx_float128/include)
target_link_libraries(test_lobatto quadmath)
add_test(NAME run_test_lobatto COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lobatto > output/test_lobatto.txt")

add_executable(test_meixner test_meixner.cpp)
target_include_directories(test_meixner PRIVATE
  ${TR29124_INCLUDE_DIRS}
  ../../wrappers)
target_link_libraries(test_meixner
  gfortran
  quadmath
  -L../../wrappers -lwrap_burkhardt)
add_test(NAME run_test_meixner COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_meixner > output/test_meixner.txt")

add_executable(test_meixner_pollaczek test_meixner_pollaczek.cpp)
target_include_directories(test_meixner_pollaczek PRIVATE )
target_link_libraries(test_meixner_pollaczek quadmath)
add_test(NAME run_test_meixner_pollaczek COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_meixner_pollaczek > output/test_meixner_pollaczek.txt")

add_executable(test_racah test_racah.cpp)
target_include_directories(test_racah PRIVATE )
target_link_libraries(test_racah quadmath)
add_test(NAME run_test_racah COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_racah > output/test_racah.txt")

add_executable(test_wilson test_wilson.cpp)
target_include_directories(test_wilson PRIVATE )
target_link_libraries(test_wilson quadmath)
add_test(NAME run_test_wilson COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_wilson > output/test_wilson.txt")

#add_executable( )
#target_include_directories( PRIVATE )
#target_link_libraries( quadmath)
