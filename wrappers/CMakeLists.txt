cmake_minimum_required (VERSION 3.2)

project(SpecfunWrappers LANGUAGES CXX C Fortran)

set(CMAKE_CXX_STANDARD 17)

set(SpecfunWrappers_VERSION_MAJOR 1)
set(SpecfunWrappers_VERSION_MINOR 0)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_C_FLAGS "-std=c11 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-fconcepts -Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_Fortran_FLAGS "-Wall -Wextra")

if (NOT MSVC)
  if (NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
endif ()

set(TR29124_LIBRARIES
  cxx_special_functions
  cxx_math_constants
  cxx_fp_utils
  cxx_complex_utils
  cxx_complex_math
  cxx_continued_fractions
  cxx_summation
  cxx_polynomial
  cxx_integration
  cxx_float128)

add_library(wrap_boost SHARED wrap_boost.cpp)
target_link_libraries(wrap_boost ${TR29124_LIBRARIES})
target_include_directories(wrap_boost PUBLIC .)


set(GSL_ROOT_DIR /usr/local)
find_package(GSL REQUIRED)

find_library(GFORTRAN_LIBRARY NAMES gfortran)

add_library(wrap_gsl SHARED
            wrap_gsl.cpp
	    ../3rdparty/gslextras/Fresnel/fresnel.c
	    ../3rdparty/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.cpp
	    ../3rdparty/gslextras/Hermite/gsl_sf_hermite.c
	   )
target_include_directories(wrap_gsl
PUBLIC
  .
PRIVATE
  ../3rdparty/gslextras)
target_link_libraries(wrap_gsl ${TR29124_LIBRARIES} gsl gslcblas)

add_library(wrap_cephes SHARED wrap_cephes.cpp)
target_include_directories(wrap_cephes PUBLIC .)
target_link_libraries(wrap_cephes ${TR29124_LIBRARIES} maths_cephes)

# This was an interesting disaster!
#	    ../cephes/cmath/sincos.c


add_library(wrap_burkhardt SHARED
            wrap_burkhardt.cpp
	    ../3rdparty/burkhardt/special_functions.f90
	    ../3rdparty/burkhardt/asa109.cpp
	    ../3rdparty/burkhardt/asa243.cpp
	    ../3rdparty/burkhardt/asa310.cpp
	    ../3rdparty/burkhardt/bernstein_polynomial.cpp
	    ../3rdparty/burkhardt/cdflib.cpp
	    ../3rdparty/burkhardt/chebyshev_polynomial.cpp
	    ../3rdparty/burkhardt/gegenbauer_polynomial.cpp
	    ../3rdparty/burkhardt/hermite_polynomial.cpp
	    ../3rdparty/burkhardt/jacobi_polynomial.cpp
	    ../3rdparty/burkhardt/kronrod.cpp
	    ../3rdparty/burkhardt/lobatto_polynomial.cpp
	    ../3rdparty/burkhardt/minmax.cpp
	    ../3rdparty/burkhardt/polpak.cpp
	    ../3rdparty/burkhardt/timestamp.cpp
	    ../3rdparty/burkhardt/toms443.cpp
	    ../3rdparty/burkhardt/toms462.cpp
	    ../3rdparty/burkhardt/toms917.cpp
	    ../3rdparty/burkhardt/quadrule.cpp
	    ../3rdparty/burkhardt/line_monte_carlo.cpp
	    ../3rdparty/burkhardt/line_ncc_rule.cpp
	    ../3rdparty/burkhardt/line_nco_rule.cpp
	    ../3rdparty/burkhardt/tetrahedron_monte_carlo.cpp
	    ../3rdparty/burkhardt/triangle_fekete_rule.cpp
	    ../3rdparty/burkhardt/triangle_felippa_rule.cpp
	    ../3rdparty/burkhardt/triangle_monte_carlo.cpp
           )
target_link_libraries(wrap_burkhardt gfortran)
target_include_directories(wrap_burkhardt
PUBLIC
  .
PRIVATE
  ../3rdparty/burkhardt)
target_link_libraries(wrap_burkhardt ${TR29124_LIBRARIES})


add_library(wrap_lerchphi SHARED
            wrap_lerch.cpp
	    ../3rdparty/lerchphi/Source/lerchphi.cpp
	   )
target_include_directories(wrap_lerchphi
PUBLIC
  .
PRIVATE
  ../3rdparty/lerchphi/Source)
target_link_libraries(wrap_lerchphi ${TR29124_LIBRARIES})


add_library(wrap_faddeeva SHARED
            wrap_faddeeva.cpp
	    ../3rdparty/Faddeeva/Faddeeva.cc
	   )
target_include_directories(wrap_faddeeva
PUBLIC
  .
PRIVATE
  ../3rdparty/Faddeeva)
target_link_libraries(wrap_faddeeva ${TR29124_LIBRARIES})


add_library(wrap_lambert SHARED
            wrap_lambert.cpp
	    ../3rdparty/LambertW/LambertW.cc
	   )
target_link_libraries(wrap_lambert ${TR29124_LIBRARIES})
target_include_directories(wrap_lambert
PUBLIC
  .
PRIVATE
  ../3rdparty/LambertW)
