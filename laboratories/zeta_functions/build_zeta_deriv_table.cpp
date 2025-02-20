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
  auto eps = emsr::epsilon(p);
  std::cerr << "eps = " << eps << '\n';

  for (int i = 10; i <= 120; ++i)
    {
      auto s = mpfr::mpreal(i, prec);
      std::cerr << "1 / s = " << 1 / s << '\n';
      int k_max = 2 + std::max(int(mpfr::pow(eps, -1 / s)), 0);
      std::cerr << "pow = " << mpfr::pow(eps, -1 / s) << '\n';
      std::cerr << "k_max = " << k_max << '\n';
      mpfr::mpreal z(0, prec);
      for (int k = k_max; k >= 2; --k)
	{
	  const auto kappa = mpfr::mpreal(k, prec);
	  z += mpfr::log(kappa) / mpfr::pow(kappa, s);
	}
      std::cout << "  " << i << "  " << std::setw(w) << z << '\n';
    }
}
