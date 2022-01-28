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
  std::size_t prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(emsr::max_digits10(p));
  auto w = 8 + std::cout.precision();

  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i, prec);
      std::cout << "  " << std::setw(w) << mpfr::zeta(x) - 1 << '\n';
    }
}
