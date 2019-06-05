cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_orthogonal_polynomials
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

enable_testing()

add_executable(debug_odd_legendre debug_odd_legendre.cpp)
target_include_directories(debug_odd_legendre PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(debug_odd_legendre quadmath)

add_executable(test_assoc_laguerre test_assoc_laguerre.cpp)
target_include_directories(test_assoc_laguerre PRIVATE )
target_link_libraries(test_assoc_laguerre quadmath)

add_executable(test_assoc_legendre test_assoc_legendre.cpp)
target_include_directories(test_assoc_legendre PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_assoc_legendre quadmath)

add_executable(test_assoc_legendre_q test_assoc_legendre_q.cpp)
target_include_directories(test_assoc_legendre_q PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_assoc_legendre_q quadmath)

add_executable(test_charlier test_charlier.cpp)
target_include_directories(test_charlier PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_charlier quadmath -L../../wrappers/debug -lwrap_burkhardt)

add_executable(test_chebyshev test_chebyshev.cpp)
target_include_directories(test_chebyshev PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_chebyshev quadmath)

add_executable(test_chebyshev_trig test_chebyshev_trig.cpp)
target_include_directories(test_chebyshev_trig PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_chebyshev_trig quadmath)

add_executable(test_chebyshev_trig_pi test_chebyshev_trig_pi.cpp)
target_include_directories(test_chebyshev_trig_pi PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_chebyshev_trig_pi quadmath)

add_executable(test_complex test_complex.cpp)
target_include_directories(test_complex PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_complex quadmath)

add_executable(test_continuous_dual_hahn test_continuous_dual_hahn.cpp)
target_include_directories(test_continuous_dual_hahn PRIVATE )
target_link_libraries(test_continuous_dual_hahn quadmath)

add_executable(test_continuous_hahn test_continuous_hahn.cpp)
target_include_directories(test_continuous_hahn PRIVATE )
target_link_libraries(test_continuous_hahn quadmath)

add_executable(test_dual_hahn test_dual_hahn.cpp)
target_include_directories(test_dual_hahn PRIVATE )
target_link_libraries(test_dual_hahn quadmath)

add_executable(test_gegenbauer test_gegenbauer.cpp)
target_include_directories(test_gegenbauer PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_gegenbauer quadmath)

add_executable(test_gegenbauer_neg_parm test_gegenbauer_neg_parm.cpp)
target_include_directories(test_gegenbauer_neg_parm PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_gegenbauer_neg_parm quadmath)

add_executable(test_gegenbauer_neg_roots test_gegenbauer_neg_roots.cpp)
target_include_directories(test_gegenbauer_neg_roots PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_gegenbauer_neg_roots quadmath)

add_executable(test_hahn test_hahn.cpp)
target_include_directories(test_hahn PRIVATE )
target_link_libraries(test_hahn quadmath)

add_executable(test_hermite test_hermite.cpp)
target_include_directories(test_hermite PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../cxx_continued_fractions/include)
target_link_libraries(test_hermite quadmath)

add_executable(test_jacobi test_jacobi.cpp)
target_include_directories(test_jacobi PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_jacobi quadmath)

add_executable(test_jacobi_neg_parm test_jacobi_neg_parm.cpp)
target_include_directories(test_jacobi_neg_parm PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_jacobi_neg_parm quadmath)

add_executable(test_jacobi_neg_roots test_jacobi_neg_roots.cpp)
target_include_directories(test_jacobi_neg_roots PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_jacobi_neg_roots quadmath)

add_executable(test_krawtchouk test_krawtchouk.cpp)
target_include_directories(test_krawtchouk PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_krawtchouk quadmath -L../../wrappers/debug -lwrap_burkhardt)

add_executable(test_laguerre test_laguerre.cpp)
target_include_directories(test_laguerre PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_laguerre quadmath)

add_executable(test_legendre test_legendre.cpp)
target_include_directories(test_legendre PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_legendre quadmath)

add_executable(test_legendre_q test_legendre_q.cpp)
target_include_directories(test_legendre_q PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_legendre_q quadmath)

add_executable(test_lobatto test_lobatto.cpp)
target_include_directories(test_lobatto PRIVATE
  .
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include)
target_link_libraries(test_lobatto quadmath)

add_executable(test_meixner test_meixner.cpp)
target_include_directories(test_meixner PRIVATE
  ../../include
  ../../cxx_fp_utils/include
  ../../polynomial/include
  ../../cxx_summation/include
  ../../quadrature/include
  ../../wrappers)
target_link_libraries(test_meixner quadmath -L../../wrappers/debug -lwrap_burkhardt)

add_executable(test_meixner_pollaczek test_meixner_pollaczek.cpp)
target_include_directories(test_meixner_pollaczek PRIVATE )
target_link_libraries(test_meixner_pollaczek quadmath)

add_executable(test_racah test_racah.cpp)
target_include_directories(test_racah PRIVATE )
target_link_libraries(test_racah quadmath)

add_executable(test_wilson test_wilson.cpp)
target_include_directories(test_wilson PRIVATE )
target_link_libraries(test_wilson quadmath)

#add_executable( )
#target_include_directories( PRIVATE )
#target_link_libraries( quadmath)
