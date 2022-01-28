/**
 *
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <emsr/sf_mod_bessel.h>

template<typename _Tnu, typename _Tp>
  void
  test_mod_bessel_asymp(_Tnu nu)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    for (int i = 1000; i <= 2000; ++i)
      {
	auto x = _Tp(i);
	auto iks = emsr::detail::cyl_bessel_ik_scaled_asymp(nu, x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << iks.I_value
		  << ' ' << std::setw(w) << iks.I_deriv
		  << ' ' << std::setw(w) << iks.K_value
		  << ' ' << std::setw(w) << iks.K_deriv
		  << ' ' << std::setw(w) << x * iks.Wronskian()
		  << '\n';
      }
  }

int
main()
{
  test_mod_bessel_asymp<double, double>(100);
}
