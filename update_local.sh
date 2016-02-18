#!/bin/bash

tool="cp -f"

${tool} math_const.h           ext

${tool} airy.tcc               bits
${tool} complex_util.h         bits
${tool} complex_util.tcc       bits
${tool} complex128.h           bits
${tool} float128.h             bits
${tool} float128.tcc           bits
${tool} numeric_limits.h       bits
${tool} sf_airy.tcc            bits
${tool} sf_bessel.tcc          bits
${tool} sf_beta.tcc            bits
${tool} sf_cardinal.tcc        bits
${tool} sf_chebyshev.tcc       bits
${tool} sf_dawson.tcc          bits
${tool} sf_ellint.tcc          bits
${tool} sf_expint.tcc          bits
${tool} sf_fresnel.tcc         bits
${tool} sf_gamma.tcc           bits
${tool} sf_gegenbauer.tcc      bits
${tool} sf_hankel.tcc          bits
${tool} sf_hankel_new.tcc      bits
${tool} sf_hermite.tcc         bits
${tool} sf_hydrogen.tcc        bits
${tool} sf_hyperg.tcc          bits
${tool} sf_hypint.tcc          bits
${tool} sf_jacobi.tcc          bits
${tool} sf_laguerre.tcc        bits
${tool} sf_legendre.tcc        bits
${tool} sf_mod_bessel.tcc      bits
${tool} sf_theta.tcc           bits
${tool} sf_trigint.tcc         bits
${tool} sf_zeta.tcc            bits  
${tool} specfun.h              bits
${tool} specfun_util.h         bits
