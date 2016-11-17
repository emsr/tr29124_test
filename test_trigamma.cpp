/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_trigamma test_trigamma.cpp -L$HOME/bin/lib64 -lquadmath
./test_trigamma > test_trigamma.txt

g++ -std=gnu++14 -Wall -Wextra -DNO_LOGBQ -I. -o test_trigamma test_trigamma.cpp -lquadmath
./test_trigamma > test_trigamma.txt

$HOME/bin/bin/g++ -std=gnu++14 -Wall -Wextra -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -DNO_LOGBQ -I. -o test_trigamma test_trigamma.cpp -lquadmath
./test_trigamma > test_trigamma.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include <bits/specfun.h>
#include "LentzContinuedFraction.tcc"

template<typename _Tp>
  _Tp
  __polygamma_series(unsigned int __n, _Tp __x)
  {
    constexpr int _S_max_iter = 1000;
    const auto _S_eps = __gnu_cxx::__epsilon(__x);
    auto __sum = _Tp{0};
    for (int __k = 0; __k < _S_max_iter; ++__k)
      {
	auto __term = std::pow(__x + _Tp(__k), _Tp(__n + 1));
	__sum += __term;
	if (std::abs(__term) < _S_eps * std::abs(__sum))
	  break;
      }
    return (__n & 1 ? +1 : -1) * __gnu_cxx::factorial<_Tp>(__n) * __sum;
  }

template<typename _Tp>
  _Tp
  __trigamma_cont_frac(_Tp __x)
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi(__x);
    const auto _S_12pi = _Tp{1} / _S_2pi / _Tp{6};

    auto __a
      = [_S_12pi](std::size_t __k, _Tp)
	{
	  if (__k == 1)
	    return _S_12pi;
	  else
	    {
	      auto __kk = _Tp(__k * __k);
	      return __kk * (__kk - _Tp{1}) / _Tp{4} / (_Tp{4} * __kk - _Tp{1});
	    }
	};
    using _AFun = decltype(__a);

    auto __b
      = [](std::size_t __k, _Tp __x)
	{ return __k == 0 ? _Tp{0} : (__k & 1 ? __x * __x : _Tp{1}); };
    using _BFun = decltype(__b);

    auto __w
      = [__a, __b](std::size_t __k, _Tp __x)
	{
	  auto __arg = _Tp{4} * __a(__k + 1, __x) / __x / __x + _Tp{1};
	  return __b(__k, __x) * (std::sqrt(__arg) - _Tp{1}) / _Tp{2};
	};
    using _WFun = decltype(__w);

    _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      G1Frac(__a, __b, __w);

    auto __g1 = G1Frac(__x);
    auto __rx = _Tp{1} / __x;

    return __rx + __rx * __rx / _Tp{2} + _S_2pi * __rx * __g1;
  }

template<typename _Tp>
  _Tp
  __tetragamma_cont_frac(_Tp __x)
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi(__x);
    const auto _S_8pi2 = _Tp{1} / _S_2pi / _S_2pi / _Tp{2};

    auto __a
      = [_S_8pi2](std::size_t __k, _Tp)
	{
	  if (__k == 1)
	    return _S_8pi2;
	  else
	    {
	      auto __j = __k & 1 ? (__k - 1) / 2 : __k / 2;
	      auto __aa = _Tp(__j) * _Tp(__j + 1) / _Tp{2} / _Tp(2 * __j + 1);
	      return (__k & 1 ? __aa * _Tp(__j + 1) : __aa * _Tp(__j));
	    }
	};
    using _AFun = decltype(__a);

    auto __b
      = [](std::size_t __k, _Tp __x)
	{ return __k == 0 ? _Tp{0} : (__k & 1 ? __x * __x : _Tp{1}); };
    using _BFun = decltype(__b);

    auto __w
      = [__a, __b](std::size_t __k, _Tp __x)
	{
	  auto __arg = _Tp{4} * __a(__k + 1, __x) / __x / __x + _Tp{1};
	  return __b(__k, __x) * (std::sqrt(__arg) - _Tp{1}) / _Tp{2};
	};
    using _WFun = decltype(__w);

    _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      G2Frac(__a, __b, __w);

    auto __g2 = G2Frac(__x);

    auto __fg = _S_2pi / __x;
    auto __xm2 = _Tp{1} / (__x * __x);

    return -__xm2 - __xm2 / __x - __fg * __fg * __g2;
  }

template<typename _Tp>
  void
  __test_trigamma(_Tp __proto)
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 200; i <= 1000; ++i)
      {
	auto x = i * (0.01);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << __trigamma_cont_frac(x)
		  << '\n';
      }

    std::cout << '\n';
    for (int i = 100; i <= 950; ++i)
      {
	auto x = i * (0.1);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << __trigamma_cont_frac(x)
		  << '\n';
      }
  }

template<typename _Tp>
  void
  __test_tetragamma(_Tp __proto)
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << '\n';
    for (int i = 200; i <= 1000; ++i)
      {
	auto x = i * (0.01);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << __tetragamma_cont_frac(x)
		  << '\n';
      }

    std::cout << '\n';
    for (int i = 100; i <= 950; ++i)
      {
	auto x = i * (0.1);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << __tetragamma_cont_frac(x)
		  << '\n';
      }
  }

int
main()
{
  __test_trigamma(1.0);
  __test_tetragamma(1.0);
}
