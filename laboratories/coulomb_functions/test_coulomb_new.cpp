/**
g++ -g -Wall -Wextra -o test_coulomb_new \
    -I../../cxx_numeric_limits/include \
    -I../../cxx_math_constants/include \
    -I../../cxx_special_functions/include \
    -I../../cxx_fp_utils/include \
    -I../../cxx_complex_utils/include \
    -I../../cxx_polynomial/include \
    test_coulomb_new.cpp

./test_coulomb_new > test_coulomb_new.txt
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <complex>

#include "sf_coulomb.h"

template<typename Tp>
  void
  test_coulomb()
  {
    for (auto lambda : {Tp{0}, Tp{0.5L}, Tp{1}})
      {
	for (auto eta : {Tp{-2}, Tp{0}, Tp{2}, Tp{10}})
	  {
	    std::cout << "\n\nlambda = " << lambda << "; eta = " << eta << '\n';
	    for (int irho = 1; irho <= 200; ++irho)
	      {
		auto rho = irho * Tp{0.1};
		auto coulomb = emsr::detail::coulomb_wave_FG(lambda, eta, rho, 0);
		std::cout << ' ' << std::setw(16) << rho
			  << ' ' << std::setw(16) << coulomb.F_value
			  << ' ' << std::setw(16) << coulomb.G_value
			  << ' ' << std::setw(16) << coulomb.F_deriv
			  << ' ' << std::setw(16) << coulomb.G_deriv
			  << ' ' << std::setw(16) << coulomb.F_exp
			  << ' ' << std::setw(16) << coulomb.G_exp
			  << '\n';
	      }
	  }
      }
  }

int
main()
{
  using Tp = double;
  test_coulomb<Tp>();
}
