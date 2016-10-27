
/*
$HOME/bin/bin/g++ -g -I. -o test_root_finding test_root_finding.cpp
./test_root_finding > test_root_finding.txt
*/

#include <ext/root_finding.h>
#include <ext/polynomial.h>
#include <iostream>
#include <iomanip>
#include <cmath>

__gnu_cxx::__root_state<double>
cos_state(double x)
{ return __gnu_cxx::__root_state<double>{std::cos(x), -std::sin(x)}; }

template<typename _Tp>
  void
  test_roots(_Tp proto = _Tp{})
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    // We need to cast to resolve the ambiguity in overloads in std::cos.
    using fun_t = _Tp (*)(_Tp);

    try
      {
	_Tp x_lower = 1.0, x_upper = 2.0, eps = 1.0e-12;
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
    catch(std::exception& err)
      {
	std::cout << "\nError: " << err.what() << '\n';
      }

    try
      {
	_Tp x_lower = 2.0, x_upper = 3.0;
	if (__gnu_cxx::__root_bracket((fun_t)&std::cos, x_lower, x_upper, 40))
	  {
	    std::cout << "cos root bracket: " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';
	  }
	else
	  std::cout << "No cos root found in " << std::setw(width) << x_lower << " <= x <= " << std::setw(width) << x_upper << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "\nError: " << err.what() << '\n';
      }

    try
      {
	_Tp x_lower = 0.0, x_upper = 30.0, eps = 1.0e-14;
	auto bracket = __gnu_cxx::__root_brackets((fun_t)&std::cos, x_lower, x_upper, 40);
	std::cout << "cos root brackets:\n";
	for (auto& br : bracket)
	  std::cout << "  " << std::setw(width) << br.first << " <= x <= " << std::setw(width) << br.second
		    << ": x = " << std::setw(width) << __gnu_cxx::__root_newton(&cos_state, br.first, br.second, eps) << '\n';
	std::cout << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "\nError: " << err.what() << '\n';
      }

    try
      {
	__gnu_cxx::_Polynomial<_Tp> P({4.0, 1.0, 2.0, 3.0});
	_Tp x_lower = -10.0, x_upper = 10.0, eps = 1.0e-14;
	auto bracket = __gnu_cxx::__root_brackets(P, x_lower, x_upper, 50 * P.degree());
	std::cout << "root brackets:\n";
	for (auto& br : bracket)
	  {
	    std::cout << "  " << std::setw(width) << br.first
	       << " <= x <= " << std::setw(width) << br.second;
	    auto r = __gnu_cxx::__root_bisect(P, br.first, br.second, eps);
	    std::cout << ": x = " << std::setw(width) << r << '\n';
	    P /= __gnu_cxx::_Polynomial<_Tp>(-r, _Tp{1});
	  }
	std::cout << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "\nError: " << err.what() << '\n';
      }
  }

int
main()
{
  test_roots<double>();
}
