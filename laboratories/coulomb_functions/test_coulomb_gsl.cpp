/**
g++ -g -Wall -Wextra -o test_coulomb_gsl -I../../cxx_numeric_limits/include \
    -I../../cxx_math_constants/include \
    -I../../cxx_special_functions/include \
    -I../../cxx_fp_utils/include \
    -I../../cxx_complex_utils/include \
    -I../../cxx_polynomial/include \
    test_coulomb_gsl.cpp -lgsl

./test_coulomb_gsl > test_coulomb_gsl.txt
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>

#include <gsl/gsl_sf_coulomb.h>

void
test_coulomb()
{
  for (auto lambda : {0.0, 0.5, 1.0})
    {
      for (auto eta : {-2.0, 0.0, 2.0, 10.0})
	{
	  std::cout << "\n\nlambda = " << lambda << "; eta = " << eta << '\n';
	  for (int irho = 1; irho <= 200; ++irho)
	    {
	      auto rho = irho * 0.1;
              const int k_lam_G = 0;
              gsl_sf_result F, Fp, G, Gp;
              double exp_F, exp_G;
              auto stat = gsl_sf_coulomb_wave_FG_e(eta, rho, lambda, k_lam_G,
                                                   &F, &Fp, &G, &Gp, &exp_F, &exp_G);
	      std::cout << ' ' << std::setw(16) << rho
			<< ' ' << std::setw(16) << F.val
			<< ' ' << std::setw(16) << G.val
			<< ' ' << std::setw(16) << Fp.val
			<< ' ' << std::setw(16) << Gp.val
			<< ' ' << std::setw(16) << exp_F
			<< ' ' << std::setw(16) << exp_G
			//<< ' ' << std::setw(2) << stat
			<< '\n';
	    }
	}
    }
}

int
main()
{
  test_coulomb();
}
