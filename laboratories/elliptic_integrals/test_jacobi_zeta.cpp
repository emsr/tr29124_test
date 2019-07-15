/**
 *
 */

#include <iostream>
#include <iomanip>
#include <bits/specfun.h>
#include <wrap_boost.h>

template<typename _Tp>
  void
  test_jacobi_zeta()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();
    const int n_phi = 40;

    _Tp _S_pi_2 = 1.5707963267948966192313216916397514L;

    for (auto k : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0})
      {
	std::cout << '\n' << '\n';
	for (int i = -n_phi; i <= +n_phi; ++i)
	  {
	    const auto phi = i * _S_pi_2 / n_phi;
	    try
	      {
		auto Z_boost = beast::jacobi_zeta(k, phi);
		auto Z_gnu = __gnu_cxx::jacobi_zeta(k, phi);
		std::cout << ' ' << std::setw(6) << k
			  << ' ' << std::setw(w) << phi
			  << ' ' << std::setw(w) << Z_gnu
			  << ' ' << std::setw(w) << Z_boost
			  << ' ' << std::setw(w) << std::abs(Z_gnu - Z_boost)
			  << '\n';
	      }
	    catch(...)
	      {
		std::cerr << "\nexception\n";
	      }
	  }
      }
  }
int
main()
{
  test_jacobi_zeta<double>();
}
