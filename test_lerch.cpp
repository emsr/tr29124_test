// $HOME/bin_specfun/bin/g++ -g -o test_lerch test_lerch.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_lerch > test_lerch.txt

// g++ -std=c++14 -g -o test_lerch test_lerch.cpp -lquadmath

// ./test_lerch > test_lerch.txt

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include "float128.h"
#include <bits/summation.h>

  //
  // A functor for a vanWijnGaarden compressor must have
  // _Tp operator()(int) that returns a term in the original defining series.
  //
  template<typename _Tp>
    class __lerch_term
    {
    public:

      using value_type = _Tp;

      __lerch_term(value_type __z, value_type __s, value_type __a)
      : _M_z{__z}, _M_s{__s}, _M_a{__a}
      { }

      value_type
      operator()(std::size_t __i) const
      {
	return std::pow(_M_z, value_type(__i))
	     / std::pow(_M_a + value_type(__i), _M_s);
      }

    private:

      value_type _M_z;
      value_type _M_s;
      value_type _M_a;
    };

  /**
   *  blows on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    __lerch_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __na = std::nearbyint(__a);
      const bool __integral = (std::abs(__a - _Tp(__na)) < _S_eps);
      if (__integral && __na <= 0)
	return _S_nan;
      //else if (std::abs(__z) >= _Tp{1})
	//throw std::domain_error("__lerch_sum: |z| > 1");
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
   */
  template<typename _Tp>
    _Tp
    __lerch_vanwijngaarden_sum(_Tp __z, _Tp __s, _Tp __a)
    {
      constexpr auto _S_nan = std::numeric_limits<_Tp>::quiet_NaN();
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      const auto __na = std::nearbyint(__a);
      const bool __integral = (std::abs(__a - _Tp(__na)) < _S_eps);
      if (__integral && __na <= _S_eps)
	return _S_nan;
      else if (std::abs(__z) >= _Tp{1})
	throw std::domain_error("__lerch_vanwijngaarden_sum: |z| > 1");
      else if (__z < _Tp{0})
	{
	  constexpr auto _S_maxit = 100000;
	  using __lerch_t = __lerch_term<_Tp>;
	  auto __lerch_fun = __lerch_t(__z, __s, __a);
	  __gnu_cxx::_VanWijngaardenSum<_Tp> __sum;
	  for (auto __k = 0; __k < _S_maxit; ++__k)
	    {
	      auto __temp = __lerch_fun(__k);
	      __sum += __temp;
	      if (std::abs(__temp / __sum) < _S_eps)
		return __sum();
	    }
	}
      else
	{
	  constexpr auto _S_maxit = 100000;
	  using __lerch_t = __lerch_term<_Tp>;
	  auto __lerch_fun = __lerch_t(__z, __s, __a);
	  __gnu_cxx::_VanWijngaardenCompressor<__lerch_t> __term(__lerch_fun);
	  __gnu_cxx::_VanWijngaardenSum<_Tp> __sum;
	  for (auto __k = 0; __k < _S_maxit; ++__k)
	    {
	      auto __temp = __term[__k];
	      __sum += __temp;
	      if (std::abs(__temp / __sum) < _S_eps)
		return __sum();
	    }
	}
    }

  /**
   *  blows on nonpositive integeral a.
   *  As usual, the binomial coefficient kills this for practical purposes.
   */
  template<typename _Tp>
    _Tp
    __lerch_double_sum(_Tp __z, _Tp __s, _Tp __a)
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
	      __gnu_cxx::_VanWijngaardenSum<_Tp> __sum(__term);
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
    __lerch(_Tp __z, _Tp __s, _Tp __a)
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
	return __lerch_sum(__z, __s, __a);
    }

int
main()
{
  using Tp = double;

  std::cout.precision(std::numeric_limits<Tp>::digits10);
  auto width = 8 + std::numeric_limits<Tp>::digits10;

  auto s = 1.0;
  auto a = 2.0;
  std::cout << std::endl;
  std::cout << " a = " << std::setw(width) << a << std::endl;
  std::cout << " s = " << std::setw(width) << s << std::endl;
  for (int iz = -99; iz <= +99; ++iz)
    {
      auto z = 0.01 * iz;
      auto lerch1 = __lerch_sum(z, s, a);
      auto lerch2 = __lerch_vanwijngaarden_sum(z, s, a);
      //auto lerch3 = __lerch_double_sum(z, s, a);
      std::cout << ' ' << std::setw(width) << z
		<< ' ' << std::setw(width) << lerch1
		<< ' ' << std::setw(width) << lerch2
		//<< ' ' << std::setw(width) << lerch3
		<< ' ' << std::setw(width) << lerch2 - lerch1
		<< std::endl;
    }

  auto z = 1.0;
  a = 1.0;
  std::cout << std::endl;
  std::cout << " z = " << std::setw(width) << z << std::endl;
  std::cout << " a = " << std::setw(width) << a << std::endl;
  for (int is = -99; is <= +99; ++is)
    {
      auto s = 0.01 * is;
      auto lerch1 = __lerch_sum(z, s, a);
      auto zeta = std::riemann_zeta(s);
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << lerch1
		<< ' ' << std::setw(width) << zeta
		<< std::endl;
    }

  std::cout << std::endl;
  for (int ia = 1; ia <= 10; ++ia)
    {
      auto a = 1.0 * ia;
      std::cout << "\n a = " << std::setw(width) << a << std::endl;
      for (int is = 0; is <= 50; ++is)
	{
	  auto s = 0.1 * is;
	  std::cout << "\n s = " << std::setw(width) << s << std::endl << std::endl;
	  for (int iz = -99; iz <= +99; ++iz)
	    {
	      auto z = 0.01 * iz;
	      auto lerch1 = __lerch_sum(z, s, a);
	      auto lerch2 = __lerch_vanwijngaarden_sum(z, s, a);
	      //auto lerch3 = __lerch_double_sum(z, s, a);
	      std::cout << ' ' << std::setw(width) << z
			<< ' ' << std::setw(width) << lerch1
			<< ' ' << std::setw(width) << lerch2
			//<< ' ' << std::setw(width) << lerch3
			<< ' ' << std::setw(width) << lerch2 - lerch1
			<< std::endl;
	    }
	}
    }

  auto lerch1 = __lerch_vanwijngaarden_sum(-0.75, Tp{1}, Tp{2});
  auto lerch2 = __lerch_vanwijngaarden_sum(-0.5, Tp{0}, Tp{1});
}
