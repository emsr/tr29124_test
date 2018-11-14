/*
$HOME/bin/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_root_finding test_root_finding.cpp -lquadmath
./test_root_finding > test_root_finding.txt

$HOME/bin/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_root_finding test_root_finding.cpp -lquadmath
./test_root_finding > test_root_finding.txt
*/

#include <ext/root_search.h>
#include <ext/polynomial.h>
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

	    auto x_steffensen = __gnu_cxx::__root_steffensen((fun_t)&std::cos, x_lower, x_upper, eps);
	    std::cout << "cos root (steffensen)    : x = " << std::setw(width) << x_steffensen << '\n';

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

/*
  template<typename _Tp>
    void
    __root_laguerre(_Polymomial<std::complex<_Tp>>& __a, std::complex<_Tp>& __x,
		    int& __its)
    {
      // Estimated fractional roundoff error.
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();

      // Number of fractional values.
      constexpr int MR = 8;
      // Fractions used to break a limit cycle.
      static const _Tp
      _S_frac[MR + 1]
      {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

      // Number of steps taken before trying a new fraction.
      constexpr int MT = 10;

      constexpr int _S_max_iter = MT * MR;

      Complex d, f;
      int __m = __a.degree();
      for (int __iter = 1;__iter <= _S_max_iter; ++__iter)
	{ // Loop over iterations up to allowed maximum.
	  __its = __iter;
	  auto __b = __a[__m];
	  auto __err = std::abs(__b);
	  std::complex<_Tp> __d{}, __f{};
	  auto __abx = std::abs(__x);
	  for (int __j = __m - 1; __j >= 0; --__j)
	    {
	      // Efficient computation of the polynomial and its first two derivatives.
	      // f stores P''(x)/2.
	      __f = __x * __f + __d;
	      __d = __x * __d + __b;
	      __b = __x * __b + __a[__j];
	      __err = __abx * __err + std::abs(__b);
	    }
	  __err *= ___S_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(__b) <= __err) // We have the root.
	    return;
	  // Use Laguerre's formula.
	  auto __g = __d / __b;
	  auto __g2 = __g * __g;
	  auto __h = __g2 - _Tp{2} * __f / __b;
	  auto __sq = std::sqrt(_Tp(__m - 1) * (_Tp(__m) * __h - __g2));
	  auto __gp = __g + __sq;
	  auto __gm = __g - __sq;
	  auto __abp = std::abs(__gp);
	  auto __abm = std::abs(__gm);
	  if (__abp < __abm)
	    __gp = __gm;
	  auto __dx = std::max(__abp, __abm) > _Tp{0}
		    ? _Tp(__m) / __gp
		    : std::polar(_Tp{1} + __abx, _Tp(__iter));
	  auto __x1 = __x - __dx;
	  if (__x == __x1)
	    return;
	  if (__iter % MT != 0)
	    __x = __x1;
	  else
	    __x -= ___S_frac[__iter / MT] * __dx;
	}
      std::__throw_runtime_error(__N("__root_laguerre: "
				     "Maximum number of iterations exceeded"));
    }

  template<typename _Tp>
    void
    qroot(_Polynomial<std::complex<_Tp>>& __p, _Tp& __b, _Tp& __c, _Tp __eps)
    {
      using _Poly = _Polynomial<std::complex<_Tp>>;
      constexpr int _S_max_iter = 20;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tiny = _Tp{100} * _S_eps;
      auto __n = __p.order();
      _Poly __q, __qq, __rem;
      for (int __iter = 0; __iter < _S_max_iter; ++__iter)
	{
	  _Poly __d(__c, __b, _Tp{1});

	  // First division: r, s.
	  divmod(__p, __d, __q, __rem);
	  auto __s = __rem[0];
	  auto __r = __rem[1];
	  // Second division: partial r, s with respect to c.
	  divmod(__q, __d, __qq, __rem);
	  auto __sc = -__rem[0];
	  auto __rc = -__rem[1];
	  auto __sb = -__c * __rc;
	  auto __rb = -__b * __rc + __sc;
	  // Solve 2x2 equation.
	  auto __dv = _Tp{1} / (__sb * __rc - __sc * __rb);
	  auto __delb = ( __r * __sc - __s * __rc) * __dv;
	  auto __delc = (-__r * __sb + __s * __rb) * __dv;
	  __b += __delb;
	  __delc = (-__r * __sb + __s * __rb) * __dv;
	  __c += __delc;
	  if ((std::abs(__delb) <= eps * std::abs(__b)
	      || std::abs(__b) < _S_tiny)
           && (std::abs(__delc) <= eps * std::abs(__c)
	      || std::abs(__c) < _S_tiny))
	    return;
	}
      std::__throw_runtime_error(__N("qroot: "
				     "Maximum number of iterations exceeded"));
    }
*/
