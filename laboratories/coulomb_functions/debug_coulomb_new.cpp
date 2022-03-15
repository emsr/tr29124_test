/*
g++ -g -Wall -Wextra -o debug_coulomb_new \
    -I../../cxx_numeric_limits/include \
    -I../../cxx_math_constants/include \
    -I../../cxx_special_functions/include \
    -I../../cxx_fp_utils/include \
    -I../../cxx_complex_utils/include \
    -I../../cxx_polynomial/include \
    debug_coulomb_new.cpp
*/

#include "sf_coulomb.h"

int
main()
{
  double lambda = 0.0;
  double eta = -2.0;
  double rho = 0.1;
  const int k_lam_G = 0;

  emsr::detail::coulomb_wave_FG(lambda, eta, rho, k_lam_G);

  rho = 2.6;
  emsr::detail::coulomb_wave_FG(lambda, eta, rho, k_lam_G);
}
