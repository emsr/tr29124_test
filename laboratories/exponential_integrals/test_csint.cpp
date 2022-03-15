/**
 *
 */

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <complex>

#include <emsr/sf_trigint.h>

int
main()
{
  std::cout.precision(8);
  std::cout.flags(std::ios::showpoint);
  const auto w = 8 + std::cout.precision();

  std::cout << ' ' << std::setw(w) << "x"
	    << ' ' << std::setw(w) << "Ci(x)"
	    << ' ' << std::setw(w) << "Si(x)"
	    << '\n';
  for (int i = 0; i <= 1000; ++i)
    {
      double x = i * 0.01;
      auto [Ci, Si] = emsr::detail::sincosint(x);
      std::cout << ' ' << std::setw(w) << x
		<< ' ' << std::setw(w) << Ci
		<< ' ' << std::setw(w) << Si
		<< '\n';
    }

  return 0;
}

