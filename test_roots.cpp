
/*
 $HOME/bin/bin/g++ -g -o test_roots test_roots.cpp

 ./test_roots > test_roots.txt

*/

#include "roots.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void
cos_state(double x, double* value, double* deriv)
{
  *value = std::cos(x);
  *deriv = -std::sin(x);
}

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  // We need to cast to resolve the ambiguity in overloads in std::cos.
  using fun_t = double (*)(double);

  {
    double x_lower = 1.0, x_upper = 2.0, eps = 1.0e-12;
    if (__gnu_cxx::__root_bracket((fun_t)&std::cos, x_lower, x_upper))
      {
	std::cout << "cos root bracket: " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';

	auto x_bisect = __gnu_cxx::__root_bisect((fun_t)&std::cos, x_lower, x_upper, eps);
	std::cout << "cos root (bisect)        : x = " << std::setw(width) << x_bisect << '\n';

	auto x_secant = __gnu_cxx::__root_secant((fun_t)&std::cos, x_lower, x_upper, eps);
	std::cout << "cos root (secant)        : x = " << std::setw(width) << x_secant << '\n';

	auto x_false_position = __gnu_cxx::__root_false_position((fun_t)&std::cos, x_lower, x_upper, eps);
	std::cout << "cos root (false position): x = " << std::setw(width) << x_false_position << '\n';

	auto x_ridder = __gnu_cxx::__root_ridder((fun_t)&std::cos, x_lower, x_upper, eps);
	std::cout << "cos root (ridder)        : x = " << std::setw(width) << x_ridder << '\n';

	auto x_brent = __gnu_cxx::__root_brent((fun_t)&std::cos, x_lower, x_upper, eps);
	std::cout << "cos root (brent)         : x = " << std::setw(width) << x_brent << '\n';

	auto x_newton = __gnu_cxx::__root_newton(&cos_state, x_lower, x_upper, eps);
	std::cout << "cos root (newton)        : x = " << std::setw(width) << x_newton << '\n';

	auto x_safe = __gnu_cxx::__root_safe(&cos_state, x_lower, x_upper, eps);
	std::cout << "cos root (safe)          : x = " << std::setw(width) << x_safe << '\n';
      }
    else
      std::cout << "No cos root found in " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';
  }

  {
    double x_lower = 2.0, x_upper = 3.0;
    if (__gnu_cxx::__root_bracket((fun_t)&std::cos, x_lower, x_upper, 40))
      {
	std::cout << "cos root bracket: " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';
      }
    else
      std::cout << "No cos root found in " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';
  }

  {
    double x_lower = 0.0, x_upper = 30.0, eps = 1.0e-14;
    auto bracket = __gnu_cxx::__root_brackets((fun_t)&std::cos, x_lower, x_upper, 40);
    std::cout << "cos root brackets:\n";
    for (auto& br : bracket)
      std::cout << "  " << std::setw(width) << br.first << " <= x <= " << std::setw(width) << br.second
		<< ": x = " << std::setw(width) << __gnu_cxx::__root_newton(&cos_state, br.first, br.second, eps) << '\n';
    std::cout << '\n';
  }
}
