/*
$HOME/bin_tr29124/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_root_finding test_root_finding.cpp -lquadmath
./test_root_finding > test_root_finding.txt

$HOME/bin/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_root_finding test_root_finding.cpp -lquadmath
./test_root_finding > test_root_finding.txt
*/

#include <ext/root_finding.h>
#include <polynomial/polynomial.h>
#include <bits/numeric_limits.h>
#include <iostream>
#include <iomanip>
#include <cmath>

template<typename _Tp>
  _Tp
  signum(_Tp x)
  {
    if (x > _Tp{0})
      return _Tp{1};
    else if (x < _Tp{0})
      return _Tp{-1};
    else
      return _Tp{0};
  }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func1_state(_Tp x)
  { return {std::pow(x, _Tp{20}) - 1, _Tp{20} * std::pow(x, _Tp{19})}; }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func2_state(_Tp x)
  {
    return {signum(x) * std::sqrt(std::abs(x)),
	    _Tp{1} / std::sqrt(std::abs(x))};
  }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func3_state(_Tp x)
  { return {x * x - _Tp{1.0e-8}, _Tp{2} * x}; }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func4_state(_Tp x)
  { return {x * std::exp(-x), std::exp(-x) - x * std::exp(-x)}; }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func5_state(_Tp x)
  {
    return {_Tp{1} / (_Tp{1} + std::exp(x)),
	    -std::exp(x) / std::pow(_Tp{1} + std::exp(x), _Tp{2})};
  }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func6_state(_Tp x)
  {
    return {std::pow(x - _Tp{1}, _Tp{7}),
	    _Tp{7} * std::pow(x - _Tp{1}, _Tp{6})};
  }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  sin_state(_Tp x)
  { return {std::sin(x), std::cos(x)}; }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  cos_state(_Tp x)
  { return {std::cos(x), -std::sin(x)}; }

template<typename _Tp>
  __gnu_cxx::__root_state<_Tp>
  func7_state(_Tp x)
  { return {-M_PI * x + M_E, -M_PI}; }

template<typename _Tp>
  void
  test_roots(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
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

	    auto x_newton = __gnu_cxx::__root_newton(cos_state<_Tp>, x_lower, x_upper, eps);
	    std::cout << "cos root (newton)        : x = " << std::setw(width) << x_newton << '\n';

	    auto x_safe = __gnu_cxx::__root_safe(cos_state<_Tp>, x_lower, x_upper, eps);
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
		    << ": x = " << std::setw(width) << __gnu_cxx::__root_newton(cos_state<_Tp>, br.first, br.second, eps) << '\n';
	std::cout << '\n';
      }
    catch(std::exception& err)
      {
	std::cout << "\nError: " << err.what() << '\n';
      }

    try
      {
	__gnu_cxx::_Polynomial<_Tp> P({-4.0, 1.0, -2.0, 3.0});
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
