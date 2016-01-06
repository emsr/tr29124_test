// $HOME/bin_specfun/bin/g++ -g -o test_lerch_transcendent test_lerch_transcendent.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_lerch_transcendent > test_lerch_transcendent.txt

// g++ -std=c++14 -g -o test_lerch_transcendent test_lerch_transcendent.cpp -lquadmath

// ./test_lerch_transcendent > test_lerch_transcendent.txt

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include "float128.h"
#include "vanWijngaardenSum.tcc"

  /**
   *  blows on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    __lerch_transcendent_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __na = std::nearbyint(__a);
      const bool __integral = (std::abs(__a - _Tp(__na)) < _S_eps);
      if (__integral && __na <= 0)
	return _S_nan;
      else if (std::abs(__z) >= _Tp{1})
	throw std::domain_error("__lerch_transcendent_sum: |z| > 1");
      else
	{
	  constexpr auto _S_maxit = 100000;
	  auto __zpow = _Tp{1};
	  auto __sum = std::pow(__a, -__s);
	  for (auto __k = 1; __k < _S_maxit; ++__k)
	    {
	      __zpow *= __z;
	      auto __term = __zpow * std::pow(__a + __k, -__s);
	      __sum += __term;
	      if (std::abs(__term / __sum) < _S_eps)
		return __sum;
	    }
	}
    }

  /**
   *  blows on nonpositive integeral a.
   *  As usual, the binomial coefficient kills this for practical purposes.
   */
  template<typename _Tp>
    _Tp
    __lerch_transcendent_double_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __na = std::nearbyint(__a);
      const bool __integral = (std::abs(__a - _Tp(__na)) < _S_eps);
      if (__integral && __na <= 0)
	return _S_nan;
      else if (__z == _Tp{1})
	return _S_nan;
      else
	{
	  constexpr auto _S_maxit = 10000;
	  auto __lerch = std::pow(__a, -__s);
	  const auto __zfrac = -__z / (_Tp{1} - __z);
	  auto __zfact = _Tp{1};
	  for (auto __n = 1; __n < _S_maxit; ++__n)
	    {
	      auto __term = std::pow(__a, -__s);
	      auto __bincoef = _Tp{1};
	      std::__detail::__vanWijngaardenSum<_Tp> __sum(__term);
	      for (auto __k = 1; __k <= __n; ++__k)
		{
		  __bincoef *= -_Tp(__n - __k + 1) / _Tp(__k);
		  __term *= __z * __bincoef * std::pow(__a + __k, -__s);
		  __sum += __term;
		}
	      __zfact *= __zfrac;
	      __lerch += __zfact * __sum();
	      if (std::abs(__zfact * __sum() / __lerch) < _S_eps)
		break;
	    }
	  __lerch /= (_Tp{1} - __z);
	  return __lerch;
	}
    }

  /**
   *  blows on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    __lerch_transcendent(_Tp __z, _Tp __s, _Tp __a)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      if (__isnan(__z) || __isnan(__s) || __isnan(__a))
	return _S_nan;

      const auto __na = std::nearbyint(__a);
      const bool __integral = (std::abs(__a - _Tp(__na)) < _S_eps);
      if (__integral && __na <= 0)
	return _S_nan;
      else
	return __lerch_transcendent_sum(__z, __s, __a);
    }

int
main()
{
  using Tp = double;

  std::cout.precision(std::numeric_limits<Tp>::digits10);
  auto width = std::numeric_limits<Tp>::digits10 + 6;

  auto s = 1.0;
  auto a = 2.0;
  std::cout << std::endl;
  std::cout << " a = " << a << std::endl;
  std::cout << " s = " << s << std::endl;
  for (int iz = -99; iz <= +99; ++iz)
    {
      auto z = 0.01 * iz;
      auto lerch1 = __lerch_transcendent_sum(z, s, a);
      auto lerch2 = __lerch_transcendent_double_sum(z, s, a);
      std::cout << ' ' << z
		<< ' ' << lerch1
		<< ' ' << lerch2
		<< std::endl;
    }

  auto z = 1.0;
  a = 1.0;
  std::cout << std::endl;
  std::cout << " z = " << z << std::endl;
  std::cout << " a = " << a << std::endl;
  for (int is = -99; is <= +99; ++is)
    {
      auto s = 0.01 * is;
      auto lerch1 = __lerch_transcendent_sum(z, s, a);
      auto zeta = std::riemann_zeta(s);
      std::cout << ' ' << s
		<< ' ' << lerch1
		<< ' ' << zeta
		<< std::endl;
    }

  std::cout << std::endl;
  for (int ia = 1; ia <= 10; ++ia)
    {
      auto a = 1.0 * ia;
      std::cout << "\n a = " << a << std::endl;
      for (int is = 0; is <= 50; ++is)
	{
	  auto s = 0.1 * is;
	  std::cout << "\n s = " << s << std::endl << std::endl;
	  for (int iz = -99; iz <= +99; ++iz)
	    {
	      auto z = 0.01 * iz;
	      auto lerch1 = __lerch_transcendent_sum(z, s, a);
	      //auto lerch2 = __lerch_transcendent_double_sum(z, s, a);
	      std::cout << ' ' << z
			<< ' ' << lerch1
			//<< ' ' << lerch2
			<< std::endl;
	    }
	}
    }
}
