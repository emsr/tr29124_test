cmake_minimum_required (VERSION 3.2)

project(SpecfunWrappers LANGUAGES CXX C Fortran)

set(SpecfunWrappers_VERSION_MAJOR 1)
set(SpecfunWrappers_VERSION_MINOR 0)

set(CMAKE_C_COMPILER $ENV{HOME}/bin/bin/gcc)
set(CMAKE_CXX_COMPILER $ENV{HOME}/bin/bin/g++)
set(CMAKE_Fortran_COMPILER $ENV{HOME}/bin/bin/gfortran)

set(CMAKE_C_FLAGS "-std=c11 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_Fortran_FLAGS "-Wall -Wextra")

add_library(wrap_boost SHARED
            wrap_boost.h wrap_boost.cpp)


set(GSL_ROOT_DIR /usr/local)
find_package(GSL REQUIRED)

add_library(wrap_gsl SHARED
            wrap_gsl.h wrap_gsl.cpp
	    gslextras/Fresnel/fresnel.c
	    gslextras/Fresnel/fresnel.h
	    gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c
	    gslextras/Hermite/gsl_sf_hermite.c)
target_link_libraries(wrap_gsl gsl)

add_library(wrap_cephes SHARED
            wrap_cephes.h wrap_cephes.cpp
	    cephes/bessel/cephes_bessel.h
	    cephes/bessel/mconf.h
	    cephes/bessel/airy.c
	    cephes/bessel/hyp2f1.c
	    cephes/bessel/hyperg.c
	    cephes/bessel/i0.c
	    cephes/bessel/i1.c
	    cephes/bessel/iv.c
	    cephes/bessel/j0.c
	    cephes/bessel/j1.c
	    cephes/bessel/jn.c
	    cephes/bessel/jv.c
	    cephes/bessel/k0.c
	    cephes/bessel/k1.c
	    cephes/bessel/kn.c
	    cephes/bessel/psi.c
	    cephes/bessel/struve.c
	    cephes/bessel/yn.c
	    cephes/cmath/cephes_cmath.h
	    cephes/cmath/mconf.h
	    cephes/cmath/acosh.c
	    cephes/cmath/asin.c
	    cephes/cmath/asinh.c
	    cephes/cmath/atan.c
	    cephes/cmath/atanh.c
	    cephes/cmath/cbrt.c
	    cephes/cmath/const.c
	    cephes/cmath/cosh.c
	    cephes/cmath/exp.c
	    cephes/cmath/exp10.c
	    cephes/cmath/exp2.c
	    cephes/cmath/fabs.c
	    cephes/cmath/floor.c
	    cephes/cmath/isnan.c
	    cephes/cmath/log.c
	    cephes/cmath/log10.c
	    cephes/cmath/log2.c
	    cephes/cmath/pow.c
	    cephes/cmath/powi.c
	    cephes/cmath/round.c
	    cephes/cmath/setprec.c
	    cephes/cmath/sin.c
	    cephes/cmath/sincos.c
	    cephes/cmath/sindg.c
	    cephes/cmath/sinh.c
	    cephes/cmath/sqrt.c
	    cephes/cmath/tan.c
	    cephes/cmath/tandg.c
	    cephes/cmath/tanh.c
	    cephes/cmath/unity.c
	    cephes/cprob/cephes_cprob.h
	    cephes/cprob/mconf.h
	    cephes/cprob/bdtr.c
	    cephes/cprob/btdtr.c
	    cephes/cprob/chdtr.c
	    cephes/cprob/drand.c
	    cephes/cprob/expx2.c
	    cephes/cprob/fdtr.c
	    cephes/cprob/gamma.c
	    cephes/cprob/gdtr.c
	    cephes/cprob/igam.c
	    cephes/cprob/igami.c
	    cephes/cprob/incbet.c
	    cephes/cprob/incbi.c
	    cephes/cprob/kolmogorov.c
	    cephes/cprob/nbdtr.c
	    cephes/cprob/ndtr.c
	    cephes/cprob/ndtri.c
	    cephes/cprob/pdtr.c
	    cephes/cprob/stdtr.c
	    cephes/ellf/cephes_ellf.h
	    cephes/ellf/mconf.h
	    cephes/ellf/ellie.c
	    cephes/ellf/ellik.c
	    cephes/ellf/ellpe.c
	    cephes/ellf/ellpj.c
	    cephes/ellf/ellpk.c
	    cephes/ldouble/cephes_ldouble.h
	    cephes/ldouble/mconf.h
	    cephes/ldouble/acoshl.c
	    cephes/ldouble/asinl.c
	    cephes/ldouble/asinhl.c
	    cephes/ldouble/atanl.c
	    cephes/ldouble/atanhl.c
	    cephes/ldouble/bdtrl.c
	    cephes/ldouble/btdtrl.c
	    cephes/ldouble/cbrtl.c
	    cephes/ldouble/coshl.c
	    cephes/ldouble/econst.c
	    cephes/ldouble/elliel.c
	    cephes/ldouble/ellikl.c
	    cephes/ldouble/ellpel.c
	    cephes/ldouble/ellpjl.c
	    cephes/ldouble/ellpkl.c
	    cephes/ldouble/expl.c
	    cephes/ldouble/expx2l.c
	    cephes/ldouble/exp10l.c
	    cephes/ldouble/exp2l.c
	    cephes/ldouble/floorl.c
	    cephes/ldouble/gammal.c
	    cephes/ldouble/gdtrl.c
	    cephes/ldouble/gelsl.c
	    cephes/ldouble/ieee.c
	    cephes/ldouble/igamil.c
	    cephes/ldouble/igaml.c
	    cephes/ldouble/incbetl.c
	    cephes/ldouble/incbil.c
	    cephes/ldouble/isnanl.c
	    cephes/ldouble/j0l.c
	    cephes/ldouble/j1l.c
	    cephes/ldouble/jnl.c
	    cephes/ldouble/logl.c
	    cephes/ldouble/log10l.c
	    cephes/ldouble/log2l.c
	    cephes/ldouble/nbdtrl.c
	    cephes/ldouble/ndtril.c
	    cephes/ldouble/ndtrl.c
	    cephes/ldouble/pdtrl.c
	    cephes/ldouble/polevll.c
	    cephes/ldouble/powl.c
	    cephes/ldouble/powil.c
	    cephes/ldouble/sinl.c
	    cephes/ldouble/sinhl.c
	    cephes/ldouble/sqrtl.c
	    cephes/ldouble/stdtrl.c
	    cephes/ldouble/tanl.c
	    cephes/ldouble/tanhl.c
	    cephes/ldouble/unityl.c
	    cephes/ldouble/ynl.c
	    cephes/misc/cephes_misc.h
	    cephes/misc/mconf.h
	    cephes/misc/beta.c
	    cephes/misc/chbevl.c
	    cephes/misc/dawsn.c
	    cephes/misc/ei.c
	    cephes/misc/expn.c
	    cephes/misc/fac.c
	    cephes/misc/fresnl.c
	    cephes/misc/mtherr.c
	    cephes/misc/planck.c
	    cephes/misc/polevl.c
	    cephes/misc/polylog.c
	    cephes/misc/rgamma.c
	    cephes/misc/shichi.c
	    cephes/misc/sici.c
	    cephes/misc/simpsn.c
	    cephes/misc/spence.c
	    cephes/misc/zeta.c
	    cephes/misc/zetac.c
	    cephes/polyn/cephes_polyn.h
	    cephes/polyn/mconf.h
	    cephes/polyn/euclid.c
	    cephes/polyn/polmisc.c
	    cephes/polyn/polrt.c
	    cephes/polyn/polyn.c
	    cephes/polyn/revers.c
	    cephes/c9x-complex/cephes_cmplx.h
	    cephes/c9x-complex/mconf.h
	    cephes/c9x-complex/cmplx.c
	    cephes/c9x-complex/cgamma.c
	    cephes/c9x-complex/cpolylog.c
           )


add_library(wrap_burkhardt SHARED
            wrap_burkhardt.h wrap_burkhardt.cpp
	    burkhardt/special_functions.f90
	    burkhardt/asa109.cpp
	    burkhardt/asa243.cpp
	    burkhardt/asa310.cpp
	    burkhardt/bernstein_polynomial.cpp
	    burkhardt/cdflib.cpp
	    burkhardt/chebyshev_polynomial.cpp
	    burkhardt/gegenbauer_polynomial.cpp
	    burkhardt/hermite_polynomial.cpp
	    burkhardt/jacobi_polynomial.cpp
	    burkhardt/kronrod.cpp
	    burkhardt/lobatto_polynomial.cpp
	    burkhardt/minmax.cpp
	    burkhardt/polpak.cpp
	    burkhardt/timestamp.cpp
	    burkhardt/toms462.cpp
	    burkhardt/quadrule.cpp
           )


add_library(wrap_faddeeva SHARED
            wrap_faddeeva.h wrap_faddeeva.cpp
	    Faddeeva/Faddeeva.h Faddeeva/Faddeeva.cpp)
