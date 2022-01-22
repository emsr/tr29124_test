/**
 *
 */

#include <mpreal.h>
#include <emsr/numeric_limits_mpreal.h>

// List of n! zeta(n+1)
int
main()
{
  std::size_t prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(emsr::max_digits10(p));
  auto w = 8 + std::cout.precision();

  mpfr::mpreal fact(1, prec);
  for (int i = 0; i <= 120; ++i)
    {
      auto x = mpfr::mpreal(i + 1, prec);
      std::cout << "  " << std::setw(w) << fact * mpfr::zeta(x) << '\n';
      fact *= (i + 1);
    }
}
