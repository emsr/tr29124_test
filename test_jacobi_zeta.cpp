/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_zeta test_jacobi_zeta.cpp -lquadmath -Lwrappers -lwrap_boost
./test_jacobi_zeta > test_jacobi_zeta.txt

*/

#include <iostream>
#include <iomanip>
#include <bits/specfun.h>
#include "wrap_boost.h"

int
main()
{
  double _S_pi_2 = 1.5707963267948966192313216916397514L;

  std::cout << '\n';
  for (auto k : {//-1.0, -0.99, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
		 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0})
    {
      std::cout << '\n';
      for (auto phi : {-_S_pi_2, -0.99 * _S_pi_2, -0.9 * _S_pi_2, -0.8 * _S_pi_2, -0.6 * _S_pi_2, -0.4 * _S_pi_2, -0.2 * _S_pi_2,
		       0.0, 0.2 * _S_pi_2, 0.4 * _S_pi_2, 0.6 * _S_pi_2, 0.8 * _S_pi_2, 0.9 * _S_pi_2, 0.99 * _S_pi_2, _S_pi_2})
	{
	  try
	    {
	      auto Lam_boost = beast::jacobi_zeta(k, phi);
	      auto Lam_gnu = __gnu_cxx::jacobi_zeta(k, phi);
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
