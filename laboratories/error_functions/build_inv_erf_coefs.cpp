/**
 *
 */

#include <mpreal.h>
#include <emsr/numeric_limits_mpreal.h>


int
main()
{
  int prec = 256;
  mpfr::mpreal p(0, prec);
  //std::cout.precision(emsr::digits10(p));
  std::cout.precision(40);
  auto w = 8 + std::cout.precision();
  std::cout << std::scientific;

  const unsigned long long n_max = 250;

  std::vector<mpfr::mpreal> a;
  mpfr::mpreal astart(1, prec);
  a.push_back(astart);
  for (unsigned long long n = 1; n < n_max; ++n)
    {
      mpfr::mpreal atemp(0, prec);
      for (unsigned long long k = 1; k <= n; ++k)
	atemp += (2 * (k - 1) + 1) * a[k - 1]
	       * (2 * (n - k) + 1) * a[n - k]
	       / (k * (2 * k - 1));
      atemp /= (2 * n + 1);
      a.push_back(atemp);
    }

  std::cout << "\n\n" << ' ' << std::setw(w) << "a_k" << '\n';
  for (auto aa : a)
    std::cout << ' ' << std::setw(w) << aa << '\n';

  std::vector<mpfr::mpreal> c;
  mpfr::mpreal cstart(1, prec);
  c.push_back(cstart);
  for (unsigned long long n = 1; n < n_max; ++n)
    {
      mpfr::mpreal ctemp(0, prec);
      for (unsigned long long k = 1; k <= n; ++k)
	ctemp += c[k - 1] * c[n - k] / k / (2 * k - 1);
      c.push_back(ctemp);
    }
  for (unsigned long long n = 1; n < n_max; ++n)
    c[n] /= (2 * n + 1);

  std::cout << "\n\n" << ' ' << std::setw(w) << "c_k" << '\n';
  for (auto cc : c)
    std::cout << ' ' << std::setw(w) << cc << '\n';
}
