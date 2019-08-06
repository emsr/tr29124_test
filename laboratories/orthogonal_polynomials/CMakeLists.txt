cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  lab_orthogonal_polynomials
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

add_custom_target(make_lab_orthogonal_polynomials_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_executable(debug_odd_legendre debug_odd_legendre.cpp)
target_link_libraries(debug_odd_legendre ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_debug_odd_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/debug_odd_legendre > output/debug_odd_legendre.txt")

add_executable(test_assoc_laguerre test_assoc_laguerre.cpp)
target_link_libraries(test_assoc_laguerre ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_assoc_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_laguerre > output/test_assoc_laguerre.txt")

add_executable(test_assoc_legendre test_assoc_legendre.cpp)
target_link_libraries(test_assoc_legendre ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_assoc_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_legendre > output/test_assoc_legendre.txt")

add_executable(test_assoc_legendre_q test_assoc_legendre_q.cpp)
target_link_libraries(test_assoc_legendre_q ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_assoc_legendre_q COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_assoc_legendre_q > output/test_assoc_legendre_q.txt")

add_executable(test_charlier test_charlier.cpp)
target_link_libraries(test_charlier
  ${TR29124_LIBRARIES}
  gfortran
  quadmath
  wrap_burkhardt)
add_test(NAME run_test_charlier COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_charlier > output/test_charlier.txt")

add_executable(test_chebyshev test_chebyshev.cpp)
target_link_libraries(test_chebyshev ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_chebyshev COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev > output/test_chebyshev.txt")

add_executable(test_chebyshev_trig test_chebyshev_trig.cpp)
target_link_libraries(test_chebyshev_trig ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_chebyshev_trig COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev_trig > output/test_chebyshev_trig.txt")

add_executable(test_chebyshev_trig_pi test_chebyshev_trig_pi.cpp)
target_link_libraries(test_chebyshev_trig_pi ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_chebyshev_trig_pi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_chebyshev_trig_pi > output/test_chebyshev_trig_pi.txt")

add_executable(test_complex test_complex.cpp)
target_link_libraries(test_complex ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_complex COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_complex > output/test_complex.txt")

add_executable(test_continuous_dual_hahn test_continuous_dual_hahn.cpp)
target_link_libraries(test_continuous_dual_hahn ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_continuous_dual_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_continuous_dual_hahn > output/test_continuous_dual_hahn.txt")

add_executable(test_continuous_hahn test_continuous_hahn.cpp)
target_link_libraries(test_continuous_hahn ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_continuous_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_continuous_hahn > output/test_continuous_hahn.txt")

add_executable(test_dual_hahn test_dual_hahn.cpp)
target_include_directories(test_dual_hahn PRIVATE)
target_link_libraries(test_dual_hahn ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_dual_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_dual_hahn > output/test_dual_hahn.txt")

add_executable(test_gegenbauer test_gegenbauer.cpp)
target_include_directories(test_gegenbauer PRIVATE
  .)
target_link_libraries(test_gegenbauer ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_gegenbauer COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer > output/test_gegenbauer.txt")

add_executable(test_gegenbauer_neg_parm test_gegenbauer_neg_parm.cpp)
target_include_directories(test_gegenbauer_neg_parm PRIVATE
  .)
target_link_libraries(test_gegenbauer_neg_parm ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_gegenbauer_neg_parm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer_neg_parm > output/test_gegenbauer_neg_parm.txt")

add_executable(test_gegenbauer_neg_roots test_gegenbauer_neg_roots.cpp)
target_include_directories(test_gegenbauer_neg_roots PRIVATE
  .)
target_link_libraries(test_gegenbauer_neg_roots ${TR29124_LIBRARIES} cxx_float128 quadmath)
# FIXME: This is broken until __float128 is back.
#add_test(NAME run_test_gegenbauer_neg_roots COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_gegenbauer_neg_roots > output/test_gegenbauer_neg_roots.txt")

add_executable(test_hahn test_hahn.cpp)
target_link_libraries(test_hahn ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_hahn COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hahn > output/test_hahn.txt")

add_executable(test_hermite test_hermite.cpp)
target_include_directories(test_hermite PRIVATE
  .)
target_link_libraries(test_hermite ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_hermite COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_hermite > output/test_hermite.txt")

add_executable(test_jacobi test_jacobi.cpp)
target_include_directories(test_jacobi PRIVATE
  .)
target_link_libraries(test_jacobi ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_jacobi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi > output/test_jacobi.txt")

add_executable(test_jacobi_neg_parm test_jacobi_neg_parm.cpp)
target_include_directories(test_jacobi_neg_parm PRIVATE
  .)
target_link_libraries(test_jacobi_neg_parm ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_jacobi_neg_parm COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_neg_parm > output/test_jacobi_neg_parm.txt")

add_executable(test_jacobi_neg_roots test_jacobi_neg_roots.cpp)
target_include_directories(test_jacobi_neg_roots PRIVATE
  .)
target_link_libraries(test_jacobi_neg_roots ${TR29124_LIBRARIES} cxx_float128 quadmath)
# FIXME: This is broken until __float128 is back.
#add_test(NAME run_test_jacobi_neg_roots COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jacobi_neg_roots > output/test_jacobi_neg_roots.txt")

add_executable(test_krawtchouk test_krawtchouk.cpp)
target_link_libraries(test_krawtchouk
  ${TR29124_LIBRARIES}
  gfortran
  quadmath
  wrap_burkhardt)
add_test(NAME run_test_krawtchouk COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_krawtchouk > output/test_krawtchouk.txt")

add_executable(test_laguerre test_laguerre.cpp)
target_link_libraries(test_laguerre ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_laguerre > output/test_laguerre.txt")

add_executable(test_legendre test_legendre.cpp)
target_link_libraries(test_legendre ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre > output/test_legendre.txt")

add_executable(test_legendre_q test_legendre_q.cpp)
target_link_libraries(test_legendre_q ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_legendre_q COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_legendre_q > output/test_legendre_q.txt")

add_executable(test_lobatto test_lobatto.cpp)
target_include_directories(test_lobatto PRIVATE
  .)
target_link_libraries(test_lobatto ${TR29124_LIBRARIES} cxx_float128 quadmath)
add_test(NAME run_test_lobatto COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_lobatto > output/test_lobatto.txt")

add_executable(test_meixner test_meixner.cpp)
target_link_libraries(test_meixner
  ${TR29124_LIBRARIES}
  gfortran
  quadmath
  wrap_burkhardt)
add_test(NAME run_test_meixner COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_meixner > output/test_meixner.txt")

add_executable(test_meixner_pollaczek test_meixner_pollaczek.cpp)
target_link_libraries(test_meixner_pollaczek ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_meixner_pollaczek COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_meixner_pollaczek > output/test_meixner_pollaczek.txt")

add_executable(test_racah test_racah.cpp)
target_link_libraries(test_racah ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_racah COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_racah > output/test_racah.txt")

add_executable(test_wilson test_wilson.cpp)
target_link_libraries(test_wilson ${TR29124_LIBRARIES} quadmath)
add_test(NAME run_test_wilson COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_wilson > output/test_wilson.txt")

#add_executable(xxx xxx.cpp)
#target_link_libraries(xxx ${TR29124_LIBRARIES} quadmath)
#add_test(NAME run_xxx COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/xxx > output/xxx.txt")

# Orthonormality tests.

add_executable(test_orthonorm_assoc_laguerre test_orthonorm_assoc_laguerre.cpp)
target_link_libraries(test_orthonorm_assoc_laguerre
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_assoc_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_assoc_laguerre > output/test_orthonorm_assoc_laguerre.txt")

add_executable(test_orthonorm_assoc_legendre test_orthonorm_assoc_legendre.cpp)
target_link_libraries(test_orthonorm_assoc_legendre
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_assoc_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_assoc_legendre > output/test_orthonorm_assoc_legendre.txt")

add_executable(test_orthonorm_sph_legendre test_orthonorm_sph_legendre.cpp)
target_link_libraries(test_orthonorm_sph_legendre
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_sph_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_sph_legendre > output/test_orthonorm_sph_legendre.txt")

add_executable(test_orthonorm_hermite test_orthonorm_hermite.cpp)
target_link_libraries(test_orthonorm_hermite
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_hermite COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_hermite > output/test_orthonorm_hermite.txt")

add_executable(test_orthonorm_laguerre test_orthonorm_laguerre.cpp)
target_link_libraries(test_orthonorm_laguerre
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_laguerre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_laguerre > output/test_orthonorm_laguerre.txt")

add_executable(test_orthonorm_legendre test_orthonorm_legendre.cpp)
target_link_libraries(test_orthonorm_legendre
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_legendre COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_legendre > output/test_orthonorm_legendre.txt")

add_executable(test_orthonorm_gegenbauer test_orthonorm_gegenbauer.cpp)
target_link_libraries(test_orthonorm_gegenbauer
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_gegenbauer COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_gegenbauer > output/test_orthonorm_gegenbauer.txt")

add_executable(test_orthonorm_jacobi test_orthonorm_jacobi.cpp)
target_link_libraries(test_orthonorm_jacobi
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_jacobi COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_jacobi > output/test_orthonorm_jacobi.txt")

add_executable(test_orthonorm_chebyshev_t test_orthonorm_chebyshev_t.cpp)
target_link_libraries(test_orthonorm_chebyshev_t
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_chebyshev_t COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_chebyshev_t > output/test_orthonorm_chebyshev_t.txt")

add_executable(test_orthonorm_chebyshev_u test_orthonorm_chebyshev_u.cpp)
target_link_libraries(test_orthonorm_chebyshev_u
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_chebyshev_u COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_chebyshev_u > output/test_orthonorm_chebyshev_u.txt")

add_executable(test_orthonorm_chebyshev_v test_orthonorm_chebyshev_v.cpp)
target_link_libraries(test_orthonorm_chebyshev_v
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_chebyshev_v COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_chebyshev_v > output/test_orthonorm_chebyshev_v.txt")

add_executable(test_orthonorm_chebyshev_w test_orthonorm_chebyshev_w.cpp)
target_link_libraries(test_orthonorm_chebyshev_w
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_chebyshev_w COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_chebyshev_w > output/test_orthonorm_chebyshev_w.txt")

add_executable(test_orthonorm_radpoly test_orthonorm_radpoly.cpp)
target_link_libraries(test_orthonorm_radpoly
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_radpoly COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_radpoly > output/test_orthonorm_radpoly.txt")

add_executable(test_orthonorm_zernike test_orthonorm_zernike.cpp)
target_link_libraries(test_orthonorm_zernike
  ${TR29124_LIBRARIES}
  quadmath)
add_test(NAME run_test_orthonorm_zernike COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_orthonorm_zernike > output/test_orthonorm_zernike.txt")
