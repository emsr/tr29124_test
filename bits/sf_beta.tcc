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

/** @file bits/sf_beta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 6, pp. 253-266
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 213-216
// (4) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.

#ifndef _GLIBCXX_BITS_SF_BETA_TCC
#define _GLIBCXX_BITS_SF_BETA_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief  Return the beta function: @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_gamma(_Tp __a, _Tp __b)
    {

      _Tp __bet;
      if (__a > __b)
	{
	  __bet = __gamma(__a) / __gamma(__a + __b);
	  __bet *= __gamma(__b);
	}
      else
	{
	  __bet = __gamma(__b) / __gamma(__a + __b);
	  __bet *= __gamma(__a);
	}

      return __bet;
    }

  /**
   * @brief  Return the beta function @f$B(a,b)@f$ using
   *	     the log gamma functions.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_lgamma(_Tp __a, _Tp __b)
    {
      _Tp __bet = __log_gamma(__a)
		+ __log_gamma(__b)
		- __log_gamma(__a + __b);
      _Tp __sign = __log_gamma_sign(__a)
		 * __log_gamma_sign(__b)
		 * __log_gamma_sign(__a + __b);

      if (__sign == _Tp{0})
        return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__bet > __gnu_cxx::__log_max<_Tp>())
        return __sign * __gnu_cxx::__infinity<_Tp>();
      else
	return __sign * std::exp(__bet);
    }


  /**
   * @brief  Return the beta function @f$B(x,y)@f$ using
   *	     the product form.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   * Here, we employ the product form:
   * @f[
   *   B(a,b) = \frac{a + b}{a b} \prod_{k=1}^{\infty}
   *            \frac{1 + (a + b) / k}{(1 + a / k) (1 + b / k)}
   *          = \frac{a + b}{ab} \prod_{k=1}^{\infty}
   *            \left[1 - \frac{ab}{(a + k)(b + k)}\right]
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta_product(_Tp __a, _Tp __b)
    {
      const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      const auto __ab = __a * __b;
      auto __bet = (__a + __b) / __ab;

      const unsigned int _S_max_iter = 100000; // Could need 1 / sqrt(_S_eps)

      const auto __apk = __a, __apb = __b;
      for (unsigned int __k = 1; __k < _S_max_iter; ++__k)
	{
	  auto __term = _Tp{1} - __ab / (++__apk) / (++__apb);
	  __bet *= __term;
	  if (std::abs(_Tp{1} - __term) < _S_eps)
	    break;
	}

      return __bet;
    }


  /**
   * @brief  Return the beta function @f$ B(a,b) @f$.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   *
   * @param __a The first argument of the beta function.
   * @param __b The second argument of the beta function.
   * @return  The beta function.
   */
  template<typename _Tp>
    _Tp
    __beta(_Tp __a, _Tp __b)
    {
      if (__isnan(__a) || __isnan(__b))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else
	return __beta_lgamma(__a, __b);
    }


  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __ibeta_cont_frac(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr unsigned int _S_itmax = 100;
      constexpr auto _S_fpmin = 1000 * __gnu_cxx::__min<_Tp>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();

      auto __apb = __a + __b;
      auto __ap1 = __a + _Tp{1};
      auto __am1 = __a - _Tp{1};
      auto __c = _Tp{1};
      auto __d = _Tp{1} - __apb * __x / __ap1;
      if (std::abs(__d) < _S_fpmin)
	__d = _S_fpmin;
      __d = _Tp{1} / __d;
      auto __h = __d;
      for (unsigned int __m = 1; __m <= _S_itmax; ++__m)
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
				 "continued fractions failed to converge"));
    }

  /**
   * Return the regularized incomplete beta function, @f$ I_x(a,b) @f$,
   * of arguments @c a, @c b, and @c x.
   *
   * The regularized incomplete beta function is defined by:
   * @f[
   *   I_x(a,b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   *
   * @param __a The first parameter
   * @param __b The second parameter
   * @param __x The argument
   */
  template<typename _Tp>
    _Tp
    __beta_inc(_Tp __a, _Tp __b, _Tp __x)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();


      if (__x < _Tp{0} || __x > _Tp{1})
	std::__throw_domain_error(__N("__beta_inc: argument out of range"));
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
	  auto __sign = __log_gamma_sign(__a + __b)
		      * __log_gamma_sign(__a) * __log_gamma_sign(__b);
	  auto __fact = __sign * std::exp(__log_gamma(__a + __b)
		      - __log_gamma(__a) - __log_gamma(__b)
		      + __a * std::log(__x) + __b * std::log(_Tp{1} - __x));

	  if (__x < (__a + _Tp{1}) / (__a + __b + _Tp{2}))
	    return __fact * __ibeta_cont_frac(__a, __b, __x) / __a;
	  else
	    return _Tp{1}
		 - __fact * __ibeta_cont_frac(__b, __a, _Tp{1} - __x) / __b;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // __GLIBCXX_BITS_SF_BETA_TCC
