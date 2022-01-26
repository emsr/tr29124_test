/**
 *
 */

#include <iostream>
#include <iomanip>

#include <emsr/specfun.h>
#include <wrap_boost.h>

int
main()
{
  double s_pi_2 = 1.5707963267948966192313216916397514L;

  std::cout << '\n';
  for (auto k : {//-1.0, -0.99, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
		 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0})
    {
      std::cout << '\n';
      for (auto phi : {-s_pi_2, -0.99 * s_pi_2, -0.9 * s_pi_2, -0.8 * s_pi_2, -0.6 * s_pi_2, -0.4 * s_pi_2, -0.2 * s_pi_2,
		       0.0, 0.2 * s_pi_2, 0.4 * s_pi_2, 0.6 * s_pi_2, 0.8 * s_pi_2, 0.9 * s_pi_2, 0.99 * s_pi_2, s_pi_2})
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
