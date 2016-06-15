#!  /bin/bash

tool="cp -f"
suffix="_specfun"

${tool} complex_ellint_rc.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rc/complex.cc 
${tool} complex_ellint_rd.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rd/complex.cc 
${tool} complex_ellint_rf.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rf/complex.cc 
${tool} complex_ellint_rg.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rg/complex.cc 
${tool} complex_ellint_rj.cc ../gcc${suffix}/libstdc++-v3/testsuite/ext/special_functions/ellint_rj/complex.cc 
