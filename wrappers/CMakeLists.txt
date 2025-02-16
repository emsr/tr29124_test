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

#find_package(Boost REQUIRED)
find_package(Boost)

if (Boost_FOUND)
  add_library(wrap_boost SHARED wrap_boost.cpp)
  target_include_directories(wrap_boost PUBLIC . ${Boost_INCLUDE_DIRS})
  message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
endif(Boost_FOUND)

set(GSL_ROOT_DIR /usr/local)
find_package(GSL REQUIRED)

find_library(GFORTRAN_LIBRARY NAMES gfortran)

add_library(wrap_gsl SHARED
            wrap_gsl.cpp
	    ../3rdparty/gslextras/Fresnel/fresnel.c
	    ../3rdparty/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.cpp)
target_include_directories(wrap_gsl
  PUBLIC .
  PRIVATE ../3rdparty/gslextras)
target_link_libraries(wrap_gsl gsl gslcblas)


add_library(wrap_gst SHARED wrap_gst.cpp)
target_include_directories(wrap_gst PUBLIC .)
target_link_libraries(wrap_gst
  gst_incomplete_gamma
  gst_scorer
  gst_toroidal_harmonic
  gst_parabolic_cylinder_2
  gst_spheroidal_harmonic)


add_library(wrap_cephes SHARED wrap_cephes.cpp)
target_include_directories(wrap_cephes PUBLIC .)
target_link_libraries(wrap_cephes maths_cephes)


add_library(wrap_burkhardt SHARED wrap_burkhardt.cpp)
target_include_directories(wrap_burkhardt PUBLIC .)
target_link_libraries(wrap_burkhardt maths_burkhardt cxx_math_constants cxx_integration)


add_library(wrap_lerchphi SHARED
            wrap_lerch.cpp
	    ../3rdparty/lerchphi/Source/lerchphi.cpp)
target_include_directories(wrap_lerchphi
  PUBLIC .
  PRIVATE ../3rdparty/lerchphi/Source)
target_link_libraries(wrap_lerchphi)


add_library(wrap_faddeeva SHARED
            wrap_faddeeva.cpp
	    ../3rdparty/Faddeeva/Faddeeva.cc)
target_include_directories(wrap_faddeeva
  PUBLIC .
  PRIVATE ../3rdparty/Faddeeva)
target_link_libraries(wrap_faddeeva cxx_math_constants)


add_library(wrap_lambert SHARED
            wrap_lambert.cpp
	    ../3rdparty/LambertW/LambertW.cc)
target_link_libraries(wrap_lambert)
target_include_directories(wrap_lambert
  PUBLIC .
  PRIVATE ../3rdparty/LambertW)
