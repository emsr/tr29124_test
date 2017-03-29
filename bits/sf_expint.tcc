// Special functions -*- C++ -*-

// Copyright (C) 2006-2017 Free Software Foundation, Inc.
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

/** @file bits/sf_expint.tcc
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
//     Ed. by Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, New-York, Section 5, pp. 228-251.
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 222-225.
//

#ifndef _GLIBCXX_BITS_SF_EXPINT_TCC
#define _GLIBCXX_BITS_SF_EXPINT_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp> _Tp __expint_E1(_Tp);

  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$ by series summation.
   * 	    This should be good for @f$ x < 1 @f$.
   *
   * The exponential integral is given by
   *  @f[
   *    E_1(x) = \int_{1}^{\infty} \frac{e^{-xt}}{t} dt
   *  @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_E1_series(_Tp __x)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      auto __term = _Tp{1};
      auto __esum = _Tp{0};
      auto __osum = _Tp{0};
      const unsigned int __max_iter = 1000;
      for (unsigned int __i = 1; __i < __max_iter; ++__i)
	{
	  __term *= - __x / __i;
	  if (std::abs(__term)
		 < _S_eps * std::min(std::abs(__esum), std::abs(__osum)))
	    break;
	  if (__term >= _Tp{0})
	    __esum += __term / __i;
	  else
	    __osum += __term / __i;
	}

      return - __esum - __osum
	     - __gnu_cxx::__const_gamma_e(__x) - std::log(__x);
    }


  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$
   * 	    by asymptotic expansion.
   *
   * The exponential integral is given by
   * @f[
   *   E_1(x) = \int_{1}^\infty \frac{e^{-xt}}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_E1_asymp(_Tp __x)
    {
      auto __term = _Tp{1};
      auto __esum = _Tp{1};
      auto __osum = _Tp{0};
      const unsigned int __max_iter = 1000;
      for (unsigned int __i = 1; __i < __max_iter; ++__i)
	{
	  auto __prev = __term;
	  __term *= - __i / __x;
	  if (std::abs(__term) > std::abs(__prev))
	    break;
	  if (__term >= _Tp{0})
	    __esum += __term;
	  else
	    __osum += __term;
	}

      return std::exp(-__x) * (__esum + __osum) / __x;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$ by series summation.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_En_series(unsigned int __n, _Tp __x)
    {
      const unsigned int _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const int __nm1 = __n - 1;
      const auto _S_gamma_E = __gnu_cxx::__const_gamma_e(__x);
      const auto __logx = std::log(__x);
      _Tp __sum = (__nm1 != 0
		? _Tp{1} / __nm1
		: -__logx - _S_gamma_E);
      _Tp __fact = _Tp{1};
      for (int __i = 1; __i <= _S_max_iter; ++__i)
	{
	  __fact *= -__x / _Tp(__i);
	  _Tp __term;
	  if ( __i != __nm1 )
	    __term = -__fact / _Tp(__i - __nm1);
	  else
	    {
	      _Tp __psi = -_S_gamma_E;
	      for (int __ii = 1; __ii <= __nm1; ++__ii)
		__psi += _Tp{1} / _Tp(__ii);
	      __term = __fact * (__psi - __logx);
	    }
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    return __sum;
	}
      std::__throw_runtime_error(__N("__expint_En_series: "
				     "series summation failed"));
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    by continued fractions.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_En_cont_frac(unsigned int __n, _Tp __x)
    {
      const unsigned int _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_fp_min = _Tp{4} * __gnu_cxx::__min(__x);
      const int __nm1 = __n - 1;
      auto __b = __x + _Tp(__n);
      auto __c = _Tp{1} / _S_fp_min;
      auto __d = _Tp{1} / __b;
      auto __h = __d;
      for ( unsigned int __i = 1; __i <= _S_max_iter; ++__i )
	{
	  auto __a = -_Tp(__i * (__nm1 + __i));
	  __b += _Tp{2};
	  __d = _Tp{1} / (__a * __d + __b);
	  if (std::abs(__d) < _S_fp_min)
	    __d = std::copysign(_S_fp_min, __d);
	  __c = __b + __a / __c;
	  if (std::abs(__c) < _S_fp_min)
	    __c = std::copysign(_S_fp_min, __c);
	  const auto __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    return __h * std::exp(-__x);
	}
      std::__throw_runtime_error(__N("__expint_En_cont_frac: "
				     "continued fraction failed"));
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$ by recursion.
   * 	    Use upward recursion for @f$ x < n @f$
   * 	    and downward recursion (Miller's algorithm) otherwise.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_En_recursion(unsigned int __n, _Tp __x)
    {
      _Tp __En;
      _Tp __E1 = __expint_E1(__x);
      if (__x < _Tp(__n))
	{
	  // Forward recursion is stable only for n < x.
	  __En = __E1;
	  for (unsigned int __j = 2; __j < __n; ++__j)
	    __En = (std::exp(-__x) - __x * __En) / _Tp(__j - 1);
	}
      else
	{
	  // Backward recursion is stable only for n >= x.
	  __En = _Tp{1};
	  /// @todo Find a principled starting number
	  /// for the @f$ E_n(x) @f$ downward recursion.
	  const int __N = __n + 20;
	  _Tp __save = _Tp{0};
	  for (int __j = __N; __j > 0; --__j)
	    {
	      __En = (std::exp(-__x) - __j * __En) / __x;
	      if (__j == __n)
		__save = __En;
	    }
	    _Tp __norm = __En / __E1;
	    __En /= __norm;
	}

      return __En;
    }

  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$ by series summation.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_Ei_series(_Tp __x)
    {
      _Tp __term = _Tp{1};
      _Tp __sum = _Tp{0};
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const unsigned int __max_iter = 1000;
      for (unsigned int __i = 1; __i < __max_iter; ++__i)
	{
	  __term *= __x / __i;
	  __sum += __term / __i;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}

      return __gnu_cxx::__const_gamma_e(__x)
	   + __sum + std::log(__x);
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$
   * 	    by asymptotic expansion.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_Ei_asymp(_Tp __x)
    {
      _Tp __term = _Tp{1};
      _Tp __sum = _Tp{1};
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const unsigned int __max_iter = 1000;
      for (unsigned int __i = 1; __i < __max_iter; ++__i)
	{
	  _Tp __prev = __term;
	  __term *= __i / __x;
	  if (std::abs(__term) >= std::abs(__prev))
	    break;
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	}

      return std::exp(__x) * __sum / __x;
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_Ei(_Tp __x)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      if (__x < _Tp{0})
	return -__expint_E1(-__x);
      else if (__x < -std::log(_S_eps))
	return __expint_Ei_series(__x);
      else
	return __expint_Ei_asymp(__x);
    }


  /**
   * @brief Return the exponential integral @f$ E_1(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   E_1(x) = \int_{1}^\infty \frac{e^{-xt}}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_E1(_Tp __x)
    {
      if (__x < _Tp{0})
	return -__expint_Ei(-__x);
      else if (__x < _Tp{1})
	return __expint_E1_series(__x);
      else if (__x < _Tp{100})
	/// @todo Find a good asymptotic switch point in @f$ E_1(x) @f$.
	return __expint_En_cont_frac(1, __x);
      else
	return __expint_E1_asymp(__x);
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    for large argument.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_En_asymp(unsigned int __n, _Tp __x)
    {
      auto __term = _Tp{1};
      auto __sum = _Tp{1};
      for (unsigned int __i = 1; __i <= __n; ++__i)
	{
	  auto __prev = __term;
	  __term *= -_Tp(__n - __i + 1) / __x;
	  if (std::abs(__term) > std::abs(__prev))
	    break;
	  __sum += __term;
	}

      return std::exp(-__x) * __sum / __x;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$
   * 	    for large order.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint_En_large_n(unsigned int __n, _Tp __x)
    {
      const auto __xpn = __x + __n;
      const auto __xpn2 = __xpn * __xpn;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      auto __term = _Tp{1};
      auto __sum = _Tp{1};
      for (unsigned int __i = 1; __i <= __n; ++__i)
	{
	  auto __prev = __term;
	  __term *= (__n - 2 * (__i - 1) * __x) / __xpn2;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	  __sum += __term;
	}

      return std::exp(-__x) * __sum / __xpn;
    }


  /**
   * @brief Return the exponential integral @f$ E_n(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   E_n(x) = \int_{1}^\infty \frac{e^{-xt}}{t^n} dt
   * @f]
   *
   * @param  __n  The order of the exponential integral function.
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint(unsigned int __n, _Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__n <= 1 && __x == _Tp{0})
	return __gnu_cxx::__infinity(__x);
      else
	{
	  if (__n == 0)
	    return std::exp(-__x) / __x;
	  else if (__n == 1)
	    return __expint_E1(__x);
	  else if (__x == _Tp{0})
	    return _Tp{1} / static_cast<_Tp>(__n - 1);
	  else if (__x < _Tp{1})
	    return __expint_En_series(__n, __x);
	  else if (__n > 50000)
	    /// @todo Study arbitrary switch to large-n @f$ E_n(x) @f$.
	    return __expint_En_large_n(__n, __x);
	  else if (__x > _Tp{100})
	    /// @todo Find a good asymptotic switch point in @f$ E_n(x) @f$.
	    return __expint_En_asymp(__n, __x);
	  else
	    return __expint_En_cont_frac(__n, __x);
	}
    }


  /**
   * @brief Return the exponential integral @f$ Ei(x) @f$.
   *
   * The exponential integral is given by
   * @f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * @f]
   *
   * @param  __x  The argument of the exponential integral function.
   * @return  The exponential integral.
   */
  template<typename _Tp>
    _Tp
    __expint(_Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else
	return __expint_Ei(__x);
    }

  /**
   * @brief Return the logarithmic integral @f$ li(x) @f$.
   *
   * The logarithmic integral is given by
   * @f[
   *   li(x) = Ei(\log(x))
   * @f]
   *
   * @param  __x  The argument of the logarithmic integral function.
   * @return  The logarithmic integral.
   */
  template<typename _Tp>
    _Tp
    __logint(const _Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (std::abs(__x) == _Tp{1})
	return __gnu_cxx::__infinity(__x);
      else
	return __expint(std::log(__x));
    }

  /**
   * @brief Return the hyperbolic cosine integral @f$ Chi(x) @f$.
   *
   * The hyperbolic cosine integral is given by
   * @f[
   *   Chi(x) = (Ei(x) - E_1(x))/ 2 = (Ei(x) + Ei(-x))/2
   * @f]
   *
   * @param  __x  The argument of the hyperbolic cosine integral function.
   * @return  The hyperbolic cosine integral.
   */
  template<typename _Tp>
    _Tp
    __coshint(const _Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__x == _Tp{0})
	return _Tp{0};
      else
	return (__expint_Ei(__x) - __expint_E1(__x)) / _Tp{2};
    }

  /**
   * @brief Return the hyperbolic sine integral @f$ Shi(x) @f$.
   *
   * The hyperbolic sine integral is given by
   * @f[
   *   Shi(x) = (Ei(x) + E_1(x))/2 = (Ei(x) - Ei(-x))/2
   * @f]
   *
   * @param  __x  The argument of the hyperbolic sine integral function.
   * @return  The hyperbolic sine integral.
   */
  template<typename _Tp>
    _Tp
    __sinhint(const _Tp __x)
    {
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else
	return (__expint_Ei(__x) + __expint_E1(__x)) / _Tp{2};
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_EXPINT_TCC
