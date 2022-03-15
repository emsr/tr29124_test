/**
 *
 */

#include <mpreal.h>
#include <vector>
#include <iostream>
#include <iomanip>

  /**
   * 
   */
std::vector<mpfr::mpreal>
gamma_reciprocal_series_coef(std::size_t n, std::size_t prec)
{
  const auto s_gamma_e = mpfr::const_euler(prec);
  auto sign = [](std::size_t i){ return (i & 1u) == 1u ? -1 : +1; };
  std::vector<mpfr::mpreal> c;
  c.push_back(mpfr::mpreal(0, prec));
  c.push_back(mpfr::mpreal(1, prec));
  for (auto j = 1u; j < n; ++j)
    {
      auto sum = mpfr::mpreal(0, prec);
      for (auto k = 1u; k < j; ++k)
	sum += sign(k) * c[k] * mpfr::zeta(mpfr::mpreal(j + 1 - k, prec));
      c.push_back((s_gamma_e * c[j] + sign(j) * sum) / j);
    }
  return c;
}

void
write_gamma_reciprocal_series_coef(std::size_t n, std::size_t prec)
{
  std::cout.precision(2 + std::numeric_limits<mpfr::mpreal>::digits10(prec));
  std::cout << std::showpoint << std::fixed;//std::scientific
  auto width = 8 + std::cout.precision();

  auto c128 = gamma_reciprocal_series_coef(n, prec);
  std::cout << '\n';
  for (auto& c : c128)
    std::cout << std::setw(width) << c << '\n';
}

int
main()
{
  write_gamma_reciprocal_series_coef(50, 128);
  write_gamma_reciprocal_series_coef(50, 256);
}
