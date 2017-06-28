/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o mpreal_gamma_recip mpreal_gamma_recip.cpp -lquadmath -lmpfr
./mpreal_gamma_recip > mpreal_gamma_recip.txt
*/

#include <mpreal.h>
#include <vector>
#include <iostream>
#include <iomanip>

  /**
   * 
   */
std::vector<mpfr::mpreal>
__gamma_reciprocal_series_coef(std::size_t __n, std::size_t prec)
{
  const auto _S_gamma_e = mpfr::const_euler(prec);
  auto __sign = [](std::size_t __i){ return (__i & 1u) == 1u ? -1 : +1; };
  std::vector<mpfr::mpreal> __c;
  __c.push_back(mpfr::mpreal(0, prec));
  __c.push_back(mpfr::mpreal(1, prec));
  for (auto __j = 1u; __j < __n; ++__j)
    {
      auto __sum = mpfr::mpreal(0, prec);
      for (auto __k = 1u; __k < __j; ++__k)
	__sum += __sign(__k) * __c[__k] * mpfr::zeta(mpfr::mpreal(__j + 1 - __k, prec));
      __c.push_back((_S_gamma_e * __c[__j] + __sign(__j) * __sum) / __j);
    }
  return __c;
}

void
write_gamma_reciprocal_series_coef(std::size_t __n, std::size_t prec)
{
  std::cout.precision(2 + std::numeric_limits<mpfr::mpreal>::digits10(prec));
  std::cout << std::showpoint << std::fixed;//std::scientific
  auto width = 8 + std::cout.precision();

  auto c128 = __gamma_reciprocal_series_coef(__n, prec);
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
