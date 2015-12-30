// Special functions -*- C++ -*-

// Copyright (C) 2006-2015 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_zeta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland
//
// References:
// (1) Handbook of Mathematical Functions,
//     Ed. by Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, New-York, Section 5, pp. 807-808.
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.

#ifndef _GLIBCXX_BITS_SF_ZETA_TCC
#define _GLIBCXX_BITS_SF_ZETA_TCC 1

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief Compute the dilogarithm function @f$ Li_2(x) @f$
   * 	    by summation for x <= 1.
   *
   * The Riemann zeta function is defined by:
   *  \f[
   * 	Li_2(x) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *  \f]
   * For |x| near 1 use the reflection formulae:
   *  \f[
   * 	Li_2(-x) + Li_2(1-x) = \frac{\pi^2}{6} - \ln(x) \ln(1-x)
   *  \f]
   *  \f[
   * 	Li_2(-x) - Li_2(1-x) - \frac{1}{2}Li_2(1-x^2) = -\frac{\pi^2}{12} - \ln(x) \ln(1-x)
   *  \f]
   * For x < 1 use the reflection formula:
   *  \f[
   * 	Li_2(1-x) - Li_2(1-\frac{1}{1-x}) - \frac{1}{2}(\ln(x))^2
   *  \f]
   */
  template<typename _Tp>
    _Tp
    __dilog(_Tp __x)
    {
      constexpr unsigned long long _S_maxit = 100000ULL;
      constexpr _Tp _S_eps = 10 * __gnu_cxx::__math_constants<_Tp>::__eps;
      constexpr _Tp _S_pipio6
	= 1.644934066848226436472415166646025189219L;
      if (__isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__x > +_Tp{1})
	std::__throw_range_error(__N("dilog: argument greater than one"));
      else if (__x < -_Tp{1})
	{
	  auto __lnfact = std::log(_Tp{1} - __x);
	  return -__dilog(_Tp{1} - _Tp{1} / (_Tp{1} - __x))
		 - _Tp{0.5L} * __lnfact * __lnfact;
	}
      else if (__x == _Tp{1})
	return _S_pipio6;
      else if (__x == -_Tp{1})
	return -_Tp{0.5L} * _S_pipio6;
      else if (__x > _Tp{0.5L})
	return _S_pipio6 - std::log(__x) * std::log(_Tp{1} - __x)
	     - __dilog(_Tp{1} - __x);
      else if (__x < -_Tp{0.5L})
	return -_Tp{0.5L} * _S_pipio6 - std::log(_Tp{1} + __x) * std::log(-__x)
	     + __dilog(_Tp{1} + __x) - __dilog(_Tp{1} - __x * __x);
      else
	{
	  _Tp __sum = 0;
	  _Tp __fact = 1;
	  for (auto __i = 1ULL; __i < _S_maxit; ++__i)
	    {
	      __fact *= __x;
	      auto __term = __fact / (__i * __i);
	      __sum += __term;
	      if (std::abs(__term) < _S_eps)
		break;
	      if (__i + 1 == _S_maxit)
		std::__throw_runtime_error("__dilog: sum failed");
	    }
	  return __sum;
	}
    }

  /**
   * @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
   * 	     by summation for s > 1.
   *
   * The Riemann zeta function is defined by:
   *  \f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *  \f]
   * For s < 1 use the reflection formula:
   *  \f[
   * 	\zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   *  \f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_sum(_Tp __s)
    {
      // A user shouldn't get to this.
      if (__s < _Tp{1})
	std::__throw_domain_error(__N("Bad argument in zeta sum."));

      constexpr unsigned int _S_max_iter = 10000;
      _Tp __zeta = _Tp{0};
      for (unsigned int __k = 1; __k < _S_max_iter; ++__k)
	{
	  _Tp __term = std::pow(static_cast<_Tp>(__k), -__s);
	  if (__term < std::numeric_limits<_Tp>::epsilon())
	    break;
	  __zeta += __term;
	}

      return __zeta;
    }


  /**
   *   @brief  Evaluate the Riemann zeta function @f$ \zeta(s) @f$
   *           by an alternate series for s > 0.
   *
   *   The Riemann zeta function is defined by:
   *    \f[
   *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *    \f]
   *   For s < 1 use the reflection formula:
   *    \f[
   *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   *    \f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_alt(_Tp __s)
    {
      constexpr auto _S_eps = __gnu_cxx::__math_constants<_Tp>::__eps;
      constexpr unsigned int _S_max_iter = 10000000;
      _Tp __sgn = _Tp{1};
      _Tp __zeta = _Tp{0};
      for (unsigned int __i = 1; __i < _S_max_iter; ++__i)
	{
	  _Tp __term = __sgn / std::pow(_Tp(__i), __s);
	  if (std::abs(__term) < _S_eps)
	    break;
	  __zeta += __term;
	  __sgn *= -_Tp{1};
	}
      __zeta /= _Tp{1} - std::pow(_Tp{2}, _Tp{1} - __s);

      return __zeta;
    }


  /**
   *   @brief  Evaluate the Riemann zeta function by series for all s != 1.
   *           Convergence is great until largish negative numbers.
   *           Then the convergence of the > 0 sum gets better.
   *
   *   The series is:
   *    \f[
   *      \zeta(s) = \frac{1}{1-2^{1-s}}
   *                 \sum_{n=0}^{\infty} \frac{1}{2^{n+1}}
   *                 \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!} (k+1)^{-s}
   *    \f]
   *   Havil 2003, p. 206.
   *
   *   The Riemann zeta function is defined by:
   *    \f[
   *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *    \f]
   *   For s < 1 use the reflection formula:
   *    \f[
   *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   *    \f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_glob(_Tp __s)
    {
      _Tp __zeta = _Tp{0};

      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      // Max e exponent before overflow.
      const auto _S_max_bincoeff = std::numeric_limits<_Tp>::max_exponent10
				 * std::log(_Tp{10}) - _Tp{1};

      // This series works until the binomial coefficient blows up
      // so use reflection.
      if (__s < _Tp{0})
	{
	  constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
	  constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
	  if (std::fmod(__s, _Tp{2}) == _Tp{0})
	    return _Tp{0};
	  else
	    {
	      _Tp __zeta = __riemann_zeta_glob(_Tp{1} - __s);
	      __zeta *= std::pow(_Tp{2} * _S_pi, __s)
		     * std::sin(_S_pi_2 * __s)
		     * std::exp(__log_gamma(_Tp{1} - __s)) / _S_pi;
	      return __zeta;
	    }
	}

      _Tp __num = _Tp{0.25L};
      const unsigned int __maxit = 10000;
      __zeta = _Tp{0.5L}; // Zeroth order contribution already calculated.
      for (unsigned int __i = 1; __i < __maxit; ++__i)
	{
	  bool __punt = false;
	  _Tp __term = _Tp{1}; // Again, the zeroth order.
	  _Tp __bincoeff = _Tp{1};
	  for (unsigned int __j = 1; __j <= __i; ++__j)
	    {
	      __bincoeff *= -_Tp(__i - __j + 1) / _Tp(__j);
	      if(std::fabs(__bincoeff) > _S_max_bincoeff )
	      {
		// This only gets hit for x << 0.
		__punt = true;
		break;
	      }
	      __term += __bincoeff * std::pow(_Tp{1 + __j}, -__s);
	    }
	  if (__punt)
	    break;
	  __term *= __num;
	  __zeta += __term;
	  if (std::abs(__term / __zeta) < _S_eps)
	    break;
	  __num *= _Tp{0.5L};
	}

      __zeta /= _Tp{1} - std::pow(_Tp{2}, _Tp{1} - __s);

      return __zeta;
    }


  /**
   *   @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
   *           using the product over prime factors.
   *    \f[
   *      \zeta(s) = \Pi_{i=1}^\infty \frac{1}{1 - p_i^{-s}}
   *    \f]
   *    where @f$ {p_i} @f$ are the prime numbers.
   *
   *   The Riemann zeta function is defined by:
   *    \f[
   *      \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *    \f]
   *   For s < 1 use the reflection formula:
   *    \f[
   *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   *    \f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_product(_Tp __s)
    {
      extern const unsigned long __prime_list[];

      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr unsigned long
      _S_num_primes = sizeof(__prime_list) / sizeof(unsigned long);

      _Tp __zeta = _Tp{1};
      for (unsigned long __i = 0; __i < _S_num_primes; ++__i)
	{
	  const _Tp __fact = _Tp{1} - std::pow(__prime_list[__i], -__s);
	  __zeta *= __fact;
	  if (_Tp{1} - __fact < _S_eps)
	    break;
	}

      __zeta = _Tp{1} / __zeta;

      return __zeta;
    }


  /**
   *   @brief  Return the Riemann zeta function @f$ \zeta(s) @f$.
   *
   *   The Riemann zeta function is defined by:
   *    \f[
   *      \zeta(s) = \sum_{k=1}^{\infty} k^{-s} for s > 1
   *                 \frac{(2\pi)^s}{pi} sin(\frac{\pi s}{2})
   *                 \Gamma (1 - s) \zeta (1 - s) for s < 1
   *    \f]
   *   For s < 1 use the reflection formula:
   *    \f[
   *      \zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   *    \f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta(_Tp __s)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      constexpr auto _S_inf = __gnu_cxx::__math_constants<_Tp>::__inf;
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__isnan(__s))
	return _S_nan;
      else if (__s == _Tp{1})
	return _S_inf;
      else if (__s < -_Tp{19})
	{
	  _Tp __zeta = __riemann_zeta_product(_Tp{1} - __s);
	  __zeta *= std::pow(_Tp{2} * _S_pi, __s)
		 * std::sin(_Tp{0.5L} * _S_pi * __s)
		 * std::exp(__log_gamma(_Tp{1} - __s))
		 / _S_pi;
	  return __zeta;
	}
      else if (__s < _Tp{20})
	{
	  // Global double sum or McLaurin?
	  bool __glob = true;
	  if (__glob)
	    return __riemann_zeta_glob(__s);
	  else
	    {
	      if (__s > _Tp{1})
		return __riemann_zeta_sum(__s);
	      else
		{
		  _Tp __zeta = std::pow(_Tp{2} * _S_pi, __s)
			     * std::sin(_Tp{0.5L} * _S_pi * __s)
			     * __gamma(_Tp{1} - __s)
			     * __riemann_zeta_sum(_Tp{1} - __s);
		  return __zeta;
		}
	    }
	}
      else
	return __riemann_zeta_product(__s);
    }


  /**
   *   @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
   *           for all s != 1 and x > -1.
   *
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta_euler_maclaurin(_Tp __s, _Tp __a)
    {
      constexpr auto _S_eps = __gnu_cxx::__math_constants<_Tp>::__eps;
      constexpr int _S_N = 10 + std::numeric_limits<_Tp>::digits10 / 2;
      constexpr int _S_jmax = 12;

      // Coefficients for Euler-Maclaurin summation:
      // B_{2j}/(2j)!
      constexpr _Tp
      _S_hzeta_c[15]
      {
	1.0000000000000000000000000000L,
	8.3333333333333333333333333333e-02L,
       -1.3888888888888888888888888889e-03L,
	3.3068783068783068783068783069e-05L,
       -8.2671957671957671957671957672e-07L,
	2.0876756987868098979210090321e-08L,
       -5.2841901386874931848476822022e-10L,
	1.3382536530684678832826980975e-11L,
       -3.3896802963225828668301953912e-13L,
	8.5860620562778445641359054504e-15L,
       -2.1748686985580618730415164239e-16L,
	5.5090028283602295152026526089e-18L,
       -1.3954464685812523340707686264e-19L,
	3.5347070396294674716932299778e-21L,
       -8.9535174270375468504026113181e-23L
      };

      const auto __pmax  = std::pow(_Tp(_S_N) + __a, -__s);
      auto __ans = __pmax * ((_S_N + __a) / (__s - _Tp{1}) + _Tp{0.5L});
      for(auto __k = 0; __k < _S_N; ++__k)
        __ans += std::pow(__k + __a, -__s);

      auto __sfact = __s;
      auto __pfact = __pmax / (_S_N + __a);
      for(auto __j = 0; __j <= _S_jmax; ++__j)
        {
	  auto __delta = _S_hzeta_c[__j + 1] * __sfact * __pfact;
	  __ans += __delta;
	  if(std::abs(__delta / __ans) < _Tp{0.5L} * _S_eps)
	    break;
	  __sfact *= (__s + _Tp(2 * __j + 1)) * (__s + _Tp(2 * __j + 2));
	  __pfact /= (_S_N + __a) * (_S_N + __a);
        }

      return __ans;
    }


  /**
   *   @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
   *           for all s != 1 and a > -1.
   *
   *   The Hurwitz zeta function is defined by:
   *   @f[
   *     \zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n + a)^s}
   *   @f]
   *   The Riemann zeta function is a special case:
   *   @f[
   *     \zeta(s) = \zeta(s,1)
   *   @f]
   *
   *   This functions uses the double sum that converges for s != 1
   *   and a > -1:
   *   @f[
   *     \zeta(s,a) = \frac{1}{s-1}
   *                \sum_{n=0}^{\infty} \frac{1}{n + 1}
   *                \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!}\frac{1}{(k+a)^s}
   *   @f]
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob(_Tp __s, _Tp __a)
    {
      constexpr _Tp _S_eps = __gnu_cxx::__math_constants<_Tp>::__eps;
      // Max before overflow?
      constexpr _Tp _S_max = __gnu_cxx::__math_constants<_Tp>::__max;

      constexpr unsigned int _S_maxit = 10000;
      // Zeroth order contribution already calculated.
      auto __zeta = _Tp{0.5L};
      for (unsigned int __n = 1; __n < _S_maxit; ++__n)
	{
	  bool __punt = false;
	  auto __term = _Tp{1}; // Again, the zeroth order.
	  auto __bincoeff = _Tp{1};
	  for (unsigned int __k = 1; __k <= __n; ++__k)
	    {
	      __bincoeff *= -_Tp(__n - __k + 1) / _Tp(__k);
	      if (std::abs(__bincoeff) > _S_max)
	      {
		__punt = true;
		break;
	      }
	      __term += __bincoeff * std::pow(_Tp(__k + __a), _Tp{1} - __s);
	    }
	  if (__punt)
	    break;
	  __term /= _Tp(__n + 1);
	  __zeta += __term;
	  if (std::abs(__term / __zeta) < _S_eps)
	    break;
	}

      __zeta /= __s - _Tp{1};

      return __zeta;
    }


  /**
   *   @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
   *           for all s != 1 and a > -1.
   *
   *   The Hurwitz zeta function is defined by:
   *   @f[
   *     \zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n+a)^s}
   *   @f]
   *   The Riemann zeta function is a special case:
   *   @f[
   *     \zeta(s) = \zeta(s,1)
   *   @f]
   */
  template<typename _Tp>
    inline _Tp
    __hurwitz_zeta(_Tp __s, _Tp __a)
    { return __hurwitz_zeta_euler_maclaurin(__s, __a); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_ZETA_TCC
