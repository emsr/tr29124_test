/**
 *
 */

#include <iostream>
#include <iomanip>

#include <emsr/special_functions.h>
#include <emsr/math_constants.h>
#include <wrap_boost.h>

int
main()
{
  double pi_2 = emsr::pi / 2.0;

  std::cout << '\n';
  for (auto k : {//-1.0, -0.99, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
		 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0})
    {
      std::cout << '\n';
      for (auto phi : {-pi_2, -0.99 * pi_2, -0.9 * pi_2, -0.8 * pi_2, -0.6 * pi_2, -0.4 * pi_2, -0.2 * pi_2,
		       0.0, 0.2 * pi_2, 0.4 * pi_2, 0.6 * pi_2, 0.8 * pi_2, 0.9 * pi_2, 0.99 * pi_2, pi_2})
	{
	  try
	    {
	      auto Lam_boost = beast::heuman_lambda(k, phi);
	      auto Lam_gnu = emsr::heuman_lambda(k, phi);
	      std::cout << ' ' << std::setw(6) << k
			<< ' ' << std::setw(12) << phi
			<< ' ' << std::setw(12) << Lam_gnu
			<< ' ' << std::setw(12) << Lam_boost
			<< ' ' << std::setw(12) << std::abs(Lam_gnu - Lam_boost)
			<< '\n';
	    }
	  catch(...)
	    {
	      std::cerr << "\nexception\n";
	    }
	}
    }
}
