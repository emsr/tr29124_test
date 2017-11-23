/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polygamma test_polygamma.cpp -L$HOME/bin/lib64 -lquadmath
./test_polygamma > test_polygamma.txt

g++ -std=gnu++14 -Wall -Wextra -I. -o test_polygamma test_polygamma.cpp -lquadmath
./test_polygamma > test_polygamma.txt

$HOME/bin/bin/g++ -std=gnu++14 -Wall -Wextra -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o test_polygamma test_polygamma.cpp -lquadmath
./test_polygamma > test_polygamma.txt
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include <bits/specfun.h>
#include <polynomial/polynomial.h>
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
  test_trigamma(_Tp __proto)
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
  test_tetragamma(_Tp __proto)
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

/**
 * Get the polygamma reflection polynomial.
 *
 * From boost (c = cos(pi*x), polynomials are even - i.e. in c*c):
 *  {-1}
 *  c * {2}
 *  { -2, -4 }
 *  c * { 16, 8 }
 *  { -16, -88, -16 }
 *  c * { 272, 416, 32 }
 *  { -272, -2880, -1824, -64 }
 *  c * { 7936, 24576, 7680, 128 }
 *  { -7936, -137216, -185856, -31616, -256 }
 *  c * { 353792, 1841152, 1304832, 128512, 512 }
 *  { -353792, -9061376, -21253376, -8728576, -518656, -1024}
 *  c * { 22368256, 175627264, 222398464, 56520704, 2084864, 2048 }
 *  // if have long long
 *  { -22368256LL, -795300864LL, -2868264960LL, -2174832640LL, -357888000LL, -8361984LL, -4096LL }
 *  c * { 1903757312LL, 21016670208LL, 41731645440LL, 20261765120LL, 2230947840LL, 33497088LL, 8192LL }
 *  { -1903757312LL, -89702612992LL, -460858269696LL, -559148810240LL, -182172651520LL, -13754155008LL, -134094848LL, -16384LL }
 *  c * { 209865342976LL, 3099269660672LL, 8885192097792LL, 7048869314560LL, 1594922762240LL, 84134068224LL, 536608768LL, 32768LL }
 *  { -209865342976LL, -12655654469632LL, -87815735738368LL, -155964390375424LL, -84842998005760LL, -13684856848384LL, -511780323328LL, -2146926592LL, -65536LL }
 *  c * { 29088885112832LL, 553753414467584LL, 2165206642589696LL, 2550316668551168LL, 985278548541440LL, 115620218667008LL, 3100738912256LL, 8588754944LL, 131072LL }
 *  { -29088885112832LL, -2184860175433728LL, -19686087844429824LL, -48165109676113920LL, -39471306959486976LL, -11124607890751488LL, -965271355195392LL, -18733264797696LL, -34357248000LL, -262144LL }
 *  c * { 4951498053124096LL, 118071834535526400LL, 603968063567560704LL, 990081991141490688LL, 584901762421358592LL, 122829335169859584LL, 7984436548730880LL, 112949304754176LL, 137433710592LL, 524288LL }
 */
template<typename _Tp>
  __gnu_cxx::_Polynomial<_Tp>
  __polygamma_poly(unsigned int __n)
  {
    __gnu_cxx::_Polynomial<_Tp> __a{_Tp{0}, _Tp{1}};
    __gnu_cxx::_Polynomial<_Tp> __b{_Tp{1}, _Tp{0}, _Tp{-1}};
    if (__n == 0)
      return __a;
    else
      {
	auto __P = __a;
	auto __Pp = __P.derivative();
	for (unsigned int __k = 0; __k < __n; ++__k)
	  {
	    __P = -(_Tp(__k + 1) * __a * __P + __b * __Pp);
	    __Pp = __P.derivative();
	  }
	return __P;
      }
  }

template<typename _Tp>
  void
  test_polygamma_poly(_Tp __proto)
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    //auto w = 8 + std::cout.precision();

    for (int n = 0; n <= 20; ++n)
      {
	auto P = __polygamma_poly<_Tp>(n);
	std::cout << P << '\n';
      }
  }


int
main()
{
  test_trigamma(1.0);

  test_tetragamma(1.0);

  test_polygamma_poly(1LL);
}
