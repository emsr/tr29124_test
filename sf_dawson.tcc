// Special functions -*- C++ -*-

// Copyright (C) 2015 Free Software Foundation, Inc.
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

/** @file bits/sf_dawson.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_SF_DAWSON_TCC
#define _GLIBCXX_SF_DAWSON_TCC 1

/// FIXME for better than float.  I'm not sure these are good even for float

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *  @brief Compute the Dawson integral using the series expansion.
   */
  template<typename _Tp>
    _Tp
    __dawson_series(_Tp __x)
    {
      auto __x2 = __x * __x;
      _Tp __sum(1);
      auto __k = 1;
      _Tp __term(1);
      while (true)
	{
	  __term *= -(_Tp{2} / _Tp(2 * __k + 1)) * __x2;
	  __sum += __term;
	  ++__k;
	  if (std::abs(__term) < std::numeric_limits<_Tp>::epsilon())
	    break;
	}
      return __x * __sum;
    }


  /**
   *  @brief Compute the Dawson integral using a sampling theorem
   *    representation.
   */
  template<typename _Tp>
    _Tp
    __dawson_const_frac(_Tp __x)
    {
      constexpr auto _S_1_sqrtpi(0.5641895835477562869480794515607726L);
      constexpr auto _S_H(0.2);
///  @todo this needs some compile-time construction!
      constexpr auto __n_max = 10;
      static _Tp __c[__n_max + 1];
      static auto __init = false;
      if (! __init)
	{
	  __init = true;
	  for (unsigned int __i = 1; __i <= __n_max; ++__i)
	    {
	      _Tp __y = _Tp(2 * __i - 1) * _S_H;
	      __c[__i] = std::exp(-__y * __y);
	    }
	}

      auto __xx = std::abs(__x);
      auto __n0 = 2 * static_cast<int>(0.5L + 0.5L * __xx / _S_H);
      auto __xp = __xx - __n0 * _S_H;
      auto __e1 = std::exp(_Tp{2} * __xp * _S_H);
      auto __e2 = __e1 * __e1;
      auto __d1 = _Tp(__n0) + _Tp{1};
      auto __d2 = __d1 - _Tp{2};
      auto __sum = _Tp{0};
      for (unsigned int __i = 1; __i <= __n_max; ++__i)
	{
	  __sum += __c[__i] * (__e1 / __d1 + _Tp{1} / (__d2 * __e1));
	  __d1 += _Tp{2};
	  __d2 -= _Tp{2};
	  __e1 *= __e2;
	}
      return std::copysign(std::exp(-__xp * __xp), __x)
	   * __sum * _S_1_sqrtpi;
    }

  /**
   *  @brief Compute the Dawson integral.
   */
  template<typename _Tp>
    _Tp
    __dawson(_Tp __x)
    {
      constexpr auto _S_nan = __gnu_cxx::__math_constants<_Tp>::__NaN;
      constexpr _Tp _S_x_min{0.2L};

      if (__isnan(__x))
	return _S_nan;
      else if (std::abs(__x) < _S_x_min)
	return __dawson_series(__x);
      else
	return __dawson_const_frac(__x);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_SF_DAWSON_TCC
