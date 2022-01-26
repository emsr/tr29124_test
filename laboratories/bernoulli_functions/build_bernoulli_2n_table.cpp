/**
 *
 */

#include <mpreal.h>
#include <emsr/numeric_limits_mpreal.h>

int
main()
{
  int prec = 128;
  mpfr::mpreal p(0, prec);
  std::cout.precision(emsr::max_digits10(p));
  auto w = 8 + std::cout.precision();
  auto s_2pi = 2 * mpfr::const_pi(prec);
  std::cout << "  2pi = " << std::setw(w) << s_2pi << '\n';

  auto fact = mpfr::mpreal(2, prec);
  for (int n = 2; n <= 200; n += 2)
    {
      if ((n / 2) % 2 == 0)
	fact *= -1;
      fact *= mpfr::mpreal(n - 1, prec) / s_2pi;
      fact *= mpfr::mpreal(n, prec) / s_2pi;

      std::cout << "  " << std::setw(w) << fact * mpfr::zeta(mpfr::mpreal(n, prec)) << '\n';
    }
}

