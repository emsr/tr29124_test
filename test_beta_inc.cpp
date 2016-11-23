/*
$HOME/bin_specfun/bin/g++ -std=c++17 -g -Wall -Wextra -Wno-psabi -I. -o test_beta_inc test_beta_inc.cpp
./test_beta_inc > test_beta_inc.txt
*/

#include <cmath>
#include <algorithm> // For clamp.
#include <limits>
#include <stdexcept>
#include <iostream>
#include <iomanip>


  //  Evaluates the continued fraction for the incomplete beta function
  //  by the modified Lentz's method
  template<typename _Tp>
    _Tp
    __ibeta_cont_frac(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr auto _S_itmax = 100;
      const auto _S_fpmin = 1000 * std::numeric_limits<_Tp>::min();
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      auto __apb = __a + __b;
      auto __ap1 = __a + _Tp{1};
      auto __am1 = __a - _Tp{1};
      auto __c = _Tp{1};
      //  First step of Lentz's method.
      auto __d = _Tp{1} - __apb * __x / __ap1;
      if (std::abs(__d) < _S_fpmin)
	__d = _S_fpmin;
      __d = _Tp{1} / __d;
      auto __h = __d;
      for (int __m = 1; __m <= _S_itmax; ++__m)
	{
	  auto __m2 = 2 * __m;

	  //  Even step of the recurrence.
	  auto __aa = _Tp(__m) * (__b - _Tp(__m)) * __x
		     / ((__am1 + _Tp(__m2)) * (__a + _Tp(__m2)));
	  __d = _Tp{1} + __aa * __d;
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = _Tp{1} + __aa / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Tp{1} / __d;
	  __h *= __d * __c;

	  //  Odd step of the recurrence.
	  __aa = -(__a + _Tp(__m)) * (__apb + _Tp(__m)) * __x
		/ ((__a + _Tp(__m2)) * (__ap1 + _Tp(__m2)));
	  __d = _Tp{1} + __aa * __d;
	  if (std::abs(__d) < _S_fpmin)
	    __d = _S_fpmin;
	  __c = _Tp{1} + __aa / __c;
	  if (std::abs(__c) < _S_fpmin)
	    __c = _S_fpmin;
	  __d = _Tp{1} / __d;
	  auto __del = __d * __c;
	  __h *= __del;

	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    return __h;
	}
      std::__throw_runtime_error(__N("__ibeta_cont_frac: "
				     "a or b too big, or _S_itmax too small"));
    }

  ///  Returns the incomplete beta function I_x(a;b).
  template<typename _Tp>
    _Tp
    __ibeta(_Tp __a, _Tp __b, _Tp __x)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);
      if (__x < _Tp{0} || __x > _Tp{1})
	std::__throw_domain_error(__N("__ibeta: argument out of range"));
      else if (__isnan(__x) || __isnan(__a) || __isnan(__b))
	return _S_NaN;
      else if (__a == _Tp{0} && __b == _Tp{0})
	return _S_NaN;
      else if (__a == _Tp{0})
	{
	  if (__x > _Tp{0})
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else if (__b == _Tp{0})
	{
	  if (__x < _Tp{1})
	    return _Tp{0};
	  else
	    return _Tp{1};
	}
      else
	{
	  _Tp __fact;
	  if (__x == _Tp{0} || __x == _Tp{1})
	    __fact = _Tp{0};
	  else
	    __fact = std::exp(std::lgamma(__a + __b)
			    - std::lgamma(__a) - std::lgamma(__b)
			    + __a * std::log(__x) + __b * std::log(_Tp{1} - __x));
	  if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2}))
	    return __fact * __ibeta_cont_frac(__a, __b, __x) / __a;
	  else
	    return _Tp{1}
		 - __fact * __ibeta_cont_frac(__b, __a, _Tp{1} - __x) / __b;
	}
    }

template<typename _Tp>
  void
  test_ibeta(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 6;

    for (int ia = 0; ia <= 10; ++ia)
      {
	auto a = _Tp(ia);
	for (int ib = 10; ib >= 0; --ib)
	  {
	    auto b = _Tp(ib);
	    std::cout << "a = " << std::setw(6) << a << '\n';
	    std::cout << "b = " << std::setw(6) << b << '\n';
	    for (int ix = 0; ix <= 100; ++ix)
	      {
		//auto x = ix * _Tp{0.01L};
		auto x = std::clamp(ix * _Tp{0.01L}, _Tp{0}, _Tp{1});
		std::cout << ' ' << std::setw(6) << x
		          << ' ' << std::setw(width) << __ibeta(a, b, x) << '\n';
	      }
	  }
      }
  }

template<typename _Tp>
  void
  stress_test_ibeta(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;

    std::cout << "a = " << std::setw(6) << _Tp{0.001L} << '\n';
    std::cout << "b = " << std::setw(6) << _Tp{20L} << '\n';
    for (int ix = 0; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.00000000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.000000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.0000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.00000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.0000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.00000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.0001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.01L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{0.001L}, _Tp{20L}, x) << '\n';
      }

    std::cout << "a = " << std::setw(6) << _Tp{20L} << '\n';
    std::cout << "b = " << std::setw(6) << _Tp{0.001L} << '\n';
    for (int ix = 0; ix < 100; ++ix)
      {
	auto x = std::clamp(ix * _Tp{0.01L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.99L} + ix * _Tp{0.0001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.9999L} + ix * _Tp{0.000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.999999L} + ix * _Tp{0.00000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.99999999L} + ix * _Tp{0.0000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.9999999999L} + ix * _Tp{0.000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.999999999999L} + ix * _Tp{0.00000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix < 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.99999999999999L} + ix * _Tp{0.0000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.9999999999999999L} + ix * _Tp{0.000000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.999999999999999999L} + ix * _Tp{0.00000000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
    for (int ix = 1; ix <= 100; ++ix)
      {
	auto x = std::clamp(_Tp{0.99999999999999999999L} + ix * _Tp{0.0000000000000000000001L}, _Tp{0}, _Tp{1});
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << __ibeta(_Tp{20L}, _Tp{0.001L}, x) << '\n';
      }
  }

int
main()
{
  test_ibeta<long double>();
  stress_test_ibeta<long double>();
}
