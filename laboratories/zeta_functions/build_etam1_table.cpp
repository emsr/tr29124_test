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

  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i, prec);
      auto pow2 = mpfr::pow(mpfr::mpreal(2, prec), 1 - x);
      std::cout << "  " << std::setw(w) << (1 - pow2) * mpfr::zeta(x) - 1 << '\n';
    }
}
