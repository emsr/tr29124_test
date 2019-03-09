cmake_minimum_required (VERSION 3.2)

project(SpecfunWrappers LANGUAGES CXX C Fortran)

set(SpecfunWrappers_VERSION_MAJOR 1)
set(SpecfunWrappers_VERSION_MINOR 0)

if (EXISTS "$ENV{HOME}/bin/bin/gfortran")
	set(CMAKE_C_COMPILER "$ENV{HOME}/bin/bin/gcc")
	set(CMAKE_CXX_COMPILER "$ENV{HOME}/bin/bin/g++")
	set(CMAKE_Fortran_COMPILER "$ENV{HOME}/bin/bin/gfortran")
elseif (EXISTS "/$ENV{MSYSTEM}/bin/bin/gfortran")
	set(CMAKE_C_COMPILER "/$ENV{MSYSTEM}/bin/gcc")
	set(CMAKE_CXX_COMPILER "/$ENV{MSYSTEM}/bin/g++")
	set(CMAKE_Fortran_COMPILER "/$ENV{MSYSTEM}/bin/gfortran")
else ()
	set(CMAKE_C_COMPILER "/usr/bin/gcc")
	set(CMAKE_CXX_COMPILER "/usr/bin/g++")
	set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
endif ()

set(CMAKE_C_FLAGS "-std=c11 -Wall -Wextra -I$ENV{HOME}/tr29124_test")
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -Wno-deprecated-declarations -Wold-style-cast -Wno-psabi -I$ENV{HOME}/tr29124_test/include -I$ENV{HOME}/tr29124_test/polynomial/include -I$ENV{HOME}/tr29124_test/quadrature/include -I$ENV{HOME}/tr29124_test/3rdparty")
set(CMAKE_Fortran_FLAGS "-Wall -Wextra -I$ENV{HOME}/tr29124_test")

add_library(wrap_boost SHARED wrap_boost.cpp)
target_include_directories(wrap_boost PRIVATE ../cxx_fp_utils/include ../cxx_summation/include)


set(GSL_ROOT_DIR /usr/local)
find_package(GSL REQUIRED)

find_library(GFORTRAN_LIBRARY NAMES gfortran)

add_library(wrap_gsl SHARED
            wrap_gsl.cpp
	    ../3rdparty/gslextras/Fresnel/fresnel.c
	    ../3rdparty/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.cpp
	    ../3rdparty/gslextras/Hermite/gsl_sf_hermite.c
	   )
target_link_libraries(wrap_gsl gsl gslcblas)
target_include_directories(wrap_gsl PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/gslextras)

add_library(wrap_cephes SHARED
            wrap_cephes.cpp
	    ../3rdparty/cephes/bessel/airy.c
	    ../3rdparty/cephes/bessel/hyp2f1.c
	    ../3rdparty/cephes/bessel/hyperg.c
	    ../3rdparty/cephes/bessel/i0.c
	    ../3rdparty/cephes/bessel/i1.c
	    ../3rdparty/cephes/bessel/iv.c
	    ../3rdparty/cephes/bessel/j0.c
	    ../3rdparty/cephes/bessel/j1.c
	    ../3rdparty/cephes/bessel/jn.c
	    ../3rdparty/cephes/bessel/jv.c
	    ../3rdparty/cephes/bessel/k0.c
	    ../3rdparty/cephes/bessel/k1.c
	    ../3rdparty/cephes/bessel/kn.c
	    ../3rdparty/cephes/bessel/psi.c
	    ../3rdparty/cephes/bessel/struve.c
	    ../3rdparty/cephes/bessel/yn.c
	    ../3rdparty/cephes/cmath/acosh.c
	    ../3rdparty/cephes/cmath/asin.c
	    ../3rdparty/cephes/cmath/asinh.c
	    ../3rdparty/cephes/cmath/atan.c
	    ../3rdparty/cephes/cmath/atanh.c
	    ../3rdparty/cephes/cmath/cbrt.c
	    ../3rdparty/cephes/cmath/const.c
	    ../3rdparty/cephes/cmath/cosh.c
	    ../3rdparty/cephes/cmath/exp.c
	    ../3rdparty/cephes/cmath/exp10.c
	    ../3rdparty/cephes/cmath/exp2.c
	    ../3rdparty/cephes/cmath/fabs.c
	    ../3rdparty/cephes/cmath/floor.c
	    ../3rdparty/cephes/cmath/isnan.c
	    ../3rdparty/cephes/cmath/log.c
	    ../3rdparty/cephes/cmath/log10.c
	    ../3rdparty/cephes/cmath/log2.c
	    ../3rdparty/cephes/cmath/pow.c
	    ../3rdparty/cephes/cmath/powi.c
	    ../3rdparty/cephes/cmath/round.c
	    ../3rdparty/cephes/cmath/setprec.c
	    ../3rdparty/cephes/cmath/sin.c
	    ../3rdparty/cephes/cmath/sincos.c
	    ../3rdparty/cephes/cmath/sindg.c
	    ../3rdparty/cephes/cmath/sinh.c
	    ../3rdparty/cephes/cmath/sqrt.c
	    ../3rdparty/cephes/cmath/tan.c
	    ../3rdparty/cephes/cmath/tandg.c
	    ../3rdparty/cephes/cmath/tanh.c
	    ../3rdparty/cephes/cmath/unity.c
	    ../3rdparty/cephes/cprob/bdtr.c
	    ../3rdparty/cephes/cprob/btdtr.c
	    ../3rdparty/cephes/cprob/chdtr.c
	    ../3rdparty/cephes/cprob/drand.c
	    ../3rdparty/cephes/cprob/expx2.c
	    ../3rdparty/cephes/cprob/fdtr.c
	    ../3rdparty/cephes/cprob/gamma.c
	    ../3rdparty/cephes/cprob/gdtr.c
	    ../3rdparty/cephes/cprob/igam.c
	    ../3rdparty/cephes/cprob/igami.c
	    ../3rdparty/cephes/cprob/incbet.c
	    ../3rdparty/cephes/cprob/incbi.c
	    ../3rdparty/cephes/cprob/kolmogorov.c
	    ../3rdparty/cephes/cprob/nbdtr.c
	    ../3rdparty/cephes/cprob/ndtr.c
	    ../3rdparty/cephes/cprob/ndtri.c
	    ../3rdparty/cephes/cprob/pdtr.c
	    ../3rdparty/cephes/cprob/stdtr.c
	    ../3rdparty/cephes/ellf/ellie.c
	    ../3rdparty/cephes/ellf/ellik.c
	    ../3rdparty/cephes/ellf/ellpe.c
	    ../3rdparty/cephes/ellf/ellpj.c
	    ../3rdparty/cephes/ellf/ellpk.c
	    ../3rdparty/cephes/ldouble/acoshl.c
	    ../3rdparty/cephes/ldouble/asinl.c
	    ../3rdparty/cephes/ldouble/asinhl.c
	    ../3rdparty/cephes/ldouble/atanl.c
	    ../3rdparty/cephes/ldouble/atanhl.c
	    ../3rdparty/cephes/ldouble/bdtrl.c
	    ../3rdparty/cephes/ldouble/btdtrl.c
	    ../3rdparty/cephes/ldouble/cbrtl.c
	    ../3rdparty/cephes/ldouble/coshl.c
	    ../3rdparty/cephes/ldouble/econst.c
	    ../3rdparty/cephes/ldouble/elliel.c
	    ../3rdparty/cephes/ldouble/ellikl.c
	    ../3rdparty/cephes/ldouble/ellpel.c
	    ../3rdparty/cephes/ldouble/ellpjl.c
	    ../3rdparty/cephes/ldouble/ellpkl.c
	    ../3rdparty/cephes/ldouble/expl.c
	    ../3rdparty/cephes/ldouble/expx2l.c
	    ../3rdparty/cephes/ldouble/exp10l.c
	    ../3rdparty/cephes/ldouble/exp2l.c
	    ../3rdparty/cephes/ldouble/floorl.c
	    ../3rdparty/cephes/ldouble/gammal.c
	    ../3rdparty/cephes/ldouble/gdtrl.c
	    ../3rdparty/cephes/ldouble/gelsl.c
	    ../3rdparty/cephes/ldouble/ieee.c
	    ../3rdparty/cephes/ldouble/igamil.c
	    ../3rdparty/cephes/ldouble/igaml.c
	    ../3rdparty/cephes/ldouble/incbetl.c
	    ../3rdparty/cephes/ldouble/incbil.c
	    ../3rdparty/cephes/ldouble/isnanl.c
	    ../3rdparty/cephes/ldouble/j0l.c
	    ../3rdparty/cephes/ldouble/j1l.c
	    ../3rdparty/cephes/ldouble/jnl.c
	    ../3rdparty/cephes/ldouble/logl.c
	    ../3rdparty/cephes/ldouble/log10l.c
	    ../3rdparty/cephes/ldouble/log2l.c
	    ../3rdparty/cephes/ldouble/nbdtrl.c
	    ../3rdparty/cephes/ldouble/ndtril.c
	    ../3rdparty/cephes/ldouble/ndtrl.c
	    ../3rdparty/cephes/ldouble/pdtrl.c
	    ../3rdparty/cephes/ldouble/polevll.c
	    ../3rdparty/cephes/ldouble/powl.c
	    ../3rdparty/cephes/ldouble/powil.c
	    ../3rdparty/cephes/ldouble/sinl.c
	    ../3rdparty/cephes/ldouble/sinhl.c
	    ../3rdparty/cephes/ldouble/sqrtl.c
	    ../3rdparty/cephes/ldouble/stdtrl.c
	    ../3rdparty/cephes/ldouble/tanl.c
	    ../3rdparty/cephes/ldouble/tanhl.c
	    ../3rdparty/cephes/ldouble/unityl.c
	    ../3rdparty/cephes/ldouble/ynl.c
	    ../3rdparty/cephes/misc/beta.c
	    ../3rdparty/cephes/misc/chbevl.c
	    ../3rdparty/cephes/misc/dawsn.c
	    ../3rdparty/cephes/misc/ei.c
	    ../3rdparty/cephes/misc/expn.c
	    ../3rdparty/cephes/misc/fac.c
	    ../3rdparty/cephes/misc/fresnl.c
	    ../3rdparty/cephes/misc/mtherr.c
	    ../3rdparty/cephes/misc/planck.c
	    ../3rdparty/cephes/misc/polevl.c
	    ../3rdparty/cephes/misc/polylog.c
	    ../3rdparty/cephes/misc/rgamma.c
	    ../3rdparty/cephes/misc/shichi.c
	    ../3rdparty/cephes/misc/sici.c
	    ../3rdparty/cephes/misc/simpsn.c
	    ../3rdparty/cephes/misc/spence.c
	    ../3rdparty/cephes/misc/zeta.c
	    ../3rdparty/cephes/misc/zetac.c
	    ../3rdparty/cephes/polyn/euclid.c
	    ../3rdparty/cephes/polyn/polmisc.c
	    ../3rdparty/cephes/polyn/polrt.c
	    ../3rdparty/cephes/polyn/polyn.c
	    ../3rdparty/cephes/polyn/revers.c
	    ../3rdparty/cephes/c9x-complex/cmplx.c
	    ../3rdparty/cephes/c9x-complex/cgamma.c
	    ../3rdparty/cephes/c9x-complex/cpolylog.c
           )
target_include_directories(wrap_cephes PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/cephes)

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
target_include_directories(wrap_burkhardt PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/burkhardt)


add_library(wrap_lerchphi SHARED
            wrap_lerch.cpp
	    ../3rdparty/lerchphi/Source/lerchphi.cpp
	   )
target_include_directories(wrap_lerchphi PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/lerchphi/Source)


add_library(wrap_faddeeva SHARED
            wrap_faddeeva.cpp
	    ../3rdparty/Faddeeva/Faddeeva.cc
	   )
target_include_directories(wrap_faddeeva PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/Faddeeva)


add_library(wrap_lambert SHARED
            wrap_lambert.cpp
	    ../3rdparty/LambertW/LambertW.cc
	   )
target_include_directories(wrap_lambert PRIVATE ../cxx_fp_utils/include ../cxx_summation/include ../3rdparty/LambertW)
