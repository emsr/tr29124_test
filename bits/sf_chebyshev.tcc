// Special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
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

/** @file bits/sf_chebyshev.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_CHEBYSHEV_TCC
#define _GLIBCXX_BITS_SF_CHEBYSHEV_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Return a Chebyshev polynomial of non-negative order @f$ n @f$
   * and real argument @f$ x @f$ by the recursion
   * @f[
   *    C_n(x) = 2xC_{n-1} - C_{n-2}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   * @param _C0 The value of the zeroth-order Chebyshev polynomial at @f$ x @f$
   * @param _C1 The value of the first-order Chebyshev polynomial at @f$ x @f$
   */
  template<typename _Tp>
    _Tp
    __chebyshev_recur(unsigned int __n, _Tp __x, _Tp _C0, _Tp _C1)
    {
      auto _Ck = _Tp{0};
      for (unsigned int __j = 1; __j < __n; ++__j)
      {
	_Ck = _Tp{2} * __x * _C1 - _C0;
	_C0 = _C1;
	_C1 = _Ck;
      }
      return _Ck;
    }

  /**
   * Return the Chebyshev polynomial of the first kind @f$ T_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the first kind is defined by:
   * @f[
   *    T_n(x) = \cos(n \theta)
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    _Tp
    __chebyshev_t(unsigned int __n, _Tp __x)
    {
      auto _T0 = _Tp{1};
      if (__n == 0)
	return _T0;

      auto _T1 = __x;
      if (__n == 1)
	return _T1;

      return __chebyshev_recur(__n, __x, _T0, _T1);
    }

  /**
   * Return the Chebyshev polynomial of the second kind @f$ U_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the second kind is defined by:
   * @f[
   *    U_n(x) = \frac{\sin \left[(n+1)\theta \right]}{\sin(\theta)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    _Tp
    __chebyshev_u(unsigned int __n, _Tp __x)
    {
      auto _U0 = _Tp{1};
      if (__n == 0)
	return _U0;

      auto _U1 = _Tp{2} * __x;
      if (__n == 1)
	return _U1;

      return __chebyshev_recur(__n, __x, _U0, _U1);
    }

  /**
   * Return the Chebyshev polynomial of the third kind @f$ V_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the third kind is defined by:
   * @f[
   *    V_n(x) = \frac{\cos \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\cos \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    _Tp
    __chebyshev_v(unsigned int __n, _Tp __x)
    {
      auto _V0 = _Tp{1};
      if (__n == 0)
	return _V0;

      auto _V1 = _Tp{2} * __x - _Tp{1};
      if (__n == 1)
	return _V1;

      return __chebyshev_recur(__n, __x, _V0, _V1);
    }

  /**
   * Return the Chebyshev polynomial of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @f$ n @f$ and real argument @f$ x @f$.
   *
   * The Chebyshev polynomial of the fourth kind is defined by:
   * @f[
   *    W_n(x) = \frac{\sin \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\sin \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam _Tp The real type of the argument
   * @param __n The non-negative integral order
   * @param __x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    _Tp
    __chebyshev_w(unsigned int __n, _Tp __x)
    {
      auto _W0 = _Tp{1};
      if (__n == 0)
	return _W0;

      auto _W1 = _Tp{2} * __x + _Tp{1};
      if (__n == 1)
	return _W1;

      return __chebyshev_recur(__n, __x, _W0, _W1);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_BITS_SF_CHEBYSHEV_TCC
