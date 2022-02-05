/**
 *
 */

#include <iostream>
#include <iomanip>

#include <mpreal.h>
#include <emsr/numeric_limits_mpreal.h>

int
main()
{
  int prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(emsr::max_digits10(p));
  auto w = 8 + std::cout.precision();
  //auto s_2pi = 2 * mpfr::const_pi(prec);
  //std::cout << "  2pi = " << std::setw(w) << s_2pi << '\n';

  for (int n = 8; n <= 128; ++n)
    {
      std::cout << "\nbits: " << std::setw(3) << n << '\n';
      for (int i = 0; i <= n; ++i)
	{
	  auto x = mpfr::pow(mpfr::mpreal(2, prec), -i);
	  std::cout << "  " << std::setw(w) << mpfr::atan(x) << '\n';
	}
    }
}
