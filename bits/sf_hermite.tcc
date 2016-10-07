// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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

/** @file bits/sf_hermite.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// Reference:
// (1) Handbook of Mathematical Functions,
//     Ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, Section 22 pp. 773-802

#ifndef _GLIBCXX_BITS_SF_HERMITE_TCC
#define _GLIBCXX_BITS_SF_HERMITE_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of order n: @f$ H_n(x) @f$ by recursion on n.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * @param __n The order of the Hermite polynomial.
   * @param __x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_recursion(unsigned int __n, _Tp __x)
    {
      // Compute H_0.
      auto __H_nm2 = _Tp{1};
      if (__n == 0)
	return __H_nm2;

      // Compute H_1.
      auto __H_nm1 = _Tp{2} * __x;
      if (__n == 1)
	return __H_nm1;

      // Compute H_n.
      _Tp __H_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __H_n = _Tp{2} * (__x * __H_nm1 - _Tp(__i - 1) * __H_nm2);
	  __H_nm2 = __H_nm1;
	  __H_nm1 = __H_n;
	}

      return __H_n;
    }

  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of large order n: @f$ H_n(x) @f$.  We assume here
   * 	    that x >= 0.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * see "Asymptotic analysis of the Hermite polynomials
   * 	  from their differential-difference equation", 
   * 	  Diego Dominici, arXiv:math/0601078v1 [math.CA] 4 Jan 2006
   * @param __n The order of the Hermite polynomial.
   * @param __x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_asymp(unsigned int __n, _Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_2pi = _S_sqrt_2
				 * __gnu_cxx::__math_constants<_Tp>::__root_pi;
      // __x >= 0 in this routine.
      const auto __xturn = std::sqrt(_Tp(2 * __n));
      if (std::abs(__x - __xturn) < _Tp{0.05L} * __xturn)
	{
	  // Transition region x ~ sqrt(2n).
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __n6th = std::pow(_Tp(__n), _Tp{1} / _Tp{6});
	  const auto __exparg = __n * std::log(__xturn) - _Tp{3} * __n_2
			      + __xturn * __x;
	  const auto __airyarg = _S_sqrt_2 * (__x - __xturn) * __n6th;
	  _Tp _Ai, _Bi, _Aip, _Bip;
	  __airy(__airyarg, _Ai, _Bi, _Aip, _Bip);
	  return _S_sqrt_2pi * __n6th * std::exp(__exparg) * _Ai;
	}
      else if (__x < __xturn)
	{
	  // Oscillatory region |x| < sqrt(2n).
	  const auto __theta = std::asin(__x / __xturn);
	  const auto __2theta = _Tp{2} * __theta;
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __exparg = __n_2 * (_Tp{2} * std::log(__xturn)
					- std::cos(__2theta));
	  const auto __arg = __theta / _Tp{2}
		+ __n_2 * (std::sin(__2theta) + __2theta - _S_pi);
	  return std::sqrt(_Tp{2} / std::cos(__theta))
	       * std::exp(__exparg) * std::cos(__arg);
	}
      else
	{
	  // Exponential region |x| > sqrt(2n).
	  const auto __sigma = std::sqrt((__x - __xturn) * (__x + __xturn));
	  const auto __exparg = _Tp{0.5L} * (__x * (__x - __sigma) - __n)
			     + __n * std::log(__sigma + __x);
	  return std::exp(__exparg)
	       * std::sqrt(_Tp{0.5L} * (_Tp{1} + __x / __sigma));
	}
    }


  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of order n: @f$ H_n(x) @f$.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * The Hermite polynomial obeys a reflection formula:
   * @f[
   *   H_n(-x) = (-1)^n H_n(x)
   * @f]
   *
   * @param __n The order of the Hermite polynomial.
   * @param __x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite(unsigned int __n, _Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__x < _Tp{0})
	return (__n % 2 == 1 ? -1 : +1) * __poly_hermite(__n, -__x);
      else if (__n > 100)
	return __poly_hermite_asymp(__n, __x);
      else
	return __poly_hermite_recursion(__n, __x);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_HERMITE_TCC
