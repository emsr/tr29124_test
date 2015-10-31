#!/bin/bash

tool="cp -f"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"

${tool} ${gcc_dir}/special_function_util.h  ./special_function_util.h
${tool} ${gcc_dir}/special_function.h	    ./special_function.h
${tool} ${gcc_dir}/bessel_function.tcc      ./bessel_function.tcc
${tool} ${gcc_dir}/beta_function.tcc	    ./beta_function.tcc
${tool} ${gcc_dir}/ell_integral.tcc	    ./ell_integral.tcc
${tool} ${gcc_dir}/exp_integral.tcc	    ./exp_integral.tcc
${tool} ${gcc_dir}/gamma.tcc		    ./gamma.tcc
${tool} ${gcc_dir}/hypergeometric.tcc	    ./hypergeometric.tcc
${tool} ${gcc_dir}/legendre_function.tcc    ./legendre_function.tcc
${tool} ${gcc_dir}/modified_bessel_func.tcc ./modified_bessel_func.tcc
${tool} ${gcc_dir}/poly_hermite.tcc	    ./poly_hermite.tcc
${tool} ${gcc_dir}/poly_laguerre.tcc	    ./poly_laguerre.tcc
${tool} ${gcc_dir}/riemann_zeta.tcc	    ./riemann_zeta.tcc
