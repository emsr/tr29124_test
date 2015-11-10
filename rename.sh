tool="svn mv"
gcc_dir="$HOME/gcc_specfun/libstdc++-v3/include/bits"
lib_dir="$HOME/gcc_specfun/libstdc++-v3"
files="include/Makefile.in
 include/bits/poly_laguerre.tcc
 include/bits/exp_integral.tcc
 include/bits/bessel_function.tcc
 include/bits/poly_hermite.tcc
 include/bits/riemann_zeta.tcc
 include/bits/beta_function.tcc
 include/bits/gamma.tcc
 include/bits/hypergeometric.tcc
 include/bits/modified_bessel_func.tcc
 include/bits/legendre_function.tcc
 include/bits/special_function_util.h
 include/bits/ell_integral.tcc
 include/Makefile.am"

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

perl -i.orig -p -e's/_GLIBCXX_BITS_BESSEL_FUNCTION_TCC/_GLIBCXX_BITS_SF_BESSEL_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_BETA_FUNCTION_TCC/_GLIBCXX_BITS_SF_BETA_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_ELL_INTEGRAL_TCC/_GLIBCXX_BITS_SF_ELLINT_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_EXP_INTEGRAL_TCC/_GLIBCXX_BITS_SF_EXPINT_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_GAMMA_TCC/_GLIBCXX_BITS_SF_GAMMA_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_HYPERGEOMETRIC_TCC/_GLIBCXX_BITS_SF_HYPERG_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_LEGENDRE_FUNCTION_TCC/_GLIBCXX_BITS_SF_LEGENDRE_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_MODIFIED_BESSEL_FUNC_TCC/_GLIBCXX_BITS_SF_MOD_BESSEL_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_POLY_HERMITE_TCC/_GLIBCXX_BITS_SF_HERMITE_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_POLY_LAGUERRE_TCC/_GLIBCXX_BITS_SF_LAGUERRE_TCC/' ${files}
perl -i.orig -p -e's/_GLIBCXX_BITS_RIEMANN_ZETA_TCC/_GLIBCXX_BITS_SF_ZETA_TCC/' ${files}

perl -i.orig -p -e's/special_function_util.h/specfun_util.h/' ${files}
perl -i.orig -p -e's/special_function.h/specfun.h/' ${files}
perl -i.orig -p -e's/bessel_function.tcc/sf_bessel.tcc/' ${files}
perl -i.orig -p -e's/beta_function.tcc/sf_beta.tcc/' ${files}
perl -i.orig -p -e's/ell_integral.tcc/sf_ellint.tcc/' ${files}
perl -i.orig -p -e's/exp_integral.tcc/sf_expint.tcc/' ${files}
perl -i.orig -p -e's/gamma.tcc/sf_gamma.tcc/' ${files}
perl -i.orig -p -e's/hypergeometric.tcc/sf_hyperg.tcc/' ${files}
perl -i.orig -p -e's/legendre_function.tcc/sf_legendre.tcc/' ${files}
perl -i.orig -p -e's/modified_bessel_func.tcc/sf_mod_bessel.tcc/' ${files}
perl -i.orig -p -e's/poly_hermite.tcc/sf_hermite.tcc/' ${files}
perl -i.orig -p -e's/poly_laguerre.tcc/sf_laguerre.tcc/' ${files}
perl -i.orig -p -e's/riemann_zeta.tcc/sf_zeta.tcc/' ${files}
