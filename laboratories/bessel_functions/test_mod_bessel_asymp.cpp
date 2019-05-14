 /*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../quadrature/include -I../../cxx_summation/include -o test_mod_bessel_asymp test_mod_bessel_asymp.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_mod_bessel_asymp

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_mod_bessel_asymp test_mod_bessel_asymp.cpp -lquadmath
PATH=$HOME/bin/lib64:$PATH ./test_mod_bessel_asymp
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <ext/cmath>

template<typename _Tnu, typename _Tp>
  void
  test_mod_bessel_asymp(_Tnu nu)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    for (int i = 1000; i <= 2000; ++i)
      {
	auto x = _Tp(i);
	auto iks = std::__detail::__cyl_bessel_ik_scaled_asymp(nu, x);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << iks.__I_value
		  << ' ' << std::setw(w) << iks.__I_deriv
		  << ' ' << std::setw(w) << iks.__K_value
		  << ' ' << std::setw(w) << iks.__K_deriv
		  << ' ' << std::setw(w) << x * iks.__Wronskian()
		  << '\n';
      }
  }

int
main()
{
  double nu = 20.0L;

  test_mod_bessel_asymp<double, double>(100);
}
