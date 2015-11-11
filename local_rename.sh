tool="mv"
gcc_dir="."

${tool} ${gcc_dir}/specfun_util.h    ${gcc_dir}/special_function_util.h 
${tool} ${gcc_dir}/specfun.h	     ${gcc_dir}/special_function.h	
${tool} ${gcc_dir}/sf_bessel.tcc     ${gcc_dir}/bessel_function.tcc	
${tool} ${gcc_dir}/sf_beta.tcc	     ${gcc_dir}/beta_function.tcc	
${tool} ${gcc_dir}/sf_ellint.tcc     ${gcc_dir}/ell_integral.tcc	
${tool} ${gcc_dir}/sf_expint.tcc     ${gcc_dir}/exp_integral.tcc	
${tool} ${gcc_dir}/sf_gamma.tcc      ${gcc_dir}/gamma.tcc		
${tool} ${gcc_dir}/sf_hyperg.tcc     ${gcc_dir}/hypergeometric.tcc	
${tool} ${gcc_dir}/sf_legendre.tcc   ${gcc_dir}/legendre_function.tcc	
${tool} ${gcc_dir}/sf_mod_bessel.tcc ${gcc_dir}/modified_bessel_func.tcc
${tool} ${gcc_dir}/sf_hermite.tcc    ${gcc_dir}/poly_hermite.tcc	
${tool} ${gcc_dir}/sf_laguerre.tcc   ${gcc_dir}/poly_laguerre.tcc	
${tool} ${gcc_dir}/sf_zeta.tcc	     ${gcc_dir}/riemann_zeta.tcc	

tool="git mv"

${tool} ${gcc_dir}/special_function_util.h  ${gcc_dir}/specfun_util.h
${tool} ${gcc_dir}/special_function.h	    ${gcc_dir}/specfun.h
${tool} ${gcc_dir}/bessel_function.tcc	    ${gcc_dir}/sf_bessel.tcc
${tool} ${gcc_dir}/beta_function.tcc	    ${gcc_dir}/sf_beta.tcc
${tool} ${gcc_dir}/ell_integral.tcc	    ${gcc_dir}/sf_ellint.tcc
${tool} ${gcc_dir}/exp_integral.tcc	    ${gcc_dir}/sf_expint.tcc
${tool} ${gcc_dir}/gamma.tcc		    ${gcc_dir}/sf_gamma.tcc
${tool} ${gcc_dir}/hypergeometric.tcc	    ${gcc_dir}/sf_hyperg.tcc
${tool} ${gcc_dir}/legendre_function.tcc    ${gcc_dir}/sf_legendre.tcc
${tool} ${gcc_dir}/modified_bessel_func.tcc ${gcc_dir}/sf_mod_bessel.tcc
${tool} ${gcc_dir}/poly_hermite.tcc	    ${gcc_dir}/sf_hermite.tcc
${tool} ${gcc_dir}/poly_laguerre.tcc	    ${gcc_dir}/sf_laguerre.tcc
${tool} ${gcc_dir}/riemann_zeta.tcc	    ${gcc_dir}/sf_zeta.tcc

perl -i.orig -p -e's/_GLIBCXX_BITS_BESSEL_FUNCTION_TCC/_GLIBCXX_BITS_SF_BESSEL_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_BETA_FUNCTION_TCC/_GLIBCXX_BITS_SF_BETA_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_ELL_INTEGRAL_TCC/_GLIBCXX_BITS_SF_ELLINT_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_EXP_INTEGRAL_TCC/_GLIBCXX_BITS_SF_EXPINT_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_GAMMA_TCC/_GLIBCXX_BITS_SF_GAMMA_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_HYPERGEOMETRIC_TCC/_GLIBCXX_BITS_SF_HYPERG_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_LEGENDRE_FUNCTION_TCC/_GLIBCXX_BITS_SF_LEGENDRE_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_MODIFIED_BESSEL_FUNC_TCC/_GLIBCXX_BITS_SF_MOD_BESSEL_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_POLY_HERMITE_TCC/_GLIBCXX_BITS_SF_HERMITE_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_POLY_LAGUERRE_TCC/_GLIBCXX_BITS_SF_LAGUERRE_TCC/' *.tcc
perl -i.orig -p -e's/_GLIBCXX_BITS_RIEMANN_ZETA_TCC/_GLIBCXX_BITS_SF_ZETA_TCC/' *.tcc

perl -i.orig -p -e's/special_function_util.h/specfun_util.h/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/special_function.h/specfun.h/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/bessel_function.tcc/sf_bessel.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/beta_function.tcc/sf_beta.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/ell_integral.tcc/sf_ellint.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/exp_integral.tcc/sf_expint.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/gamma.tcc/sf_gamma.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/hypergeometric.tcc/sf_hyperg.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/legendre_function.tcc/sf_legendre.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/modified_bessel_func.tcc/sf_mod_bessel.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/poly_hermite.tcc/sf_hermite.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/poly_laguerre.tcc/sf_laguerre.tcc/'  diff_std.sh fetch_std.sh
perl -i.orig -p -e's/riemann_zeta.tcc/sf_zeta.tcc/'  diff_std.sh fetch_std.sh


