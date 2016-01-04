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

/** @file bits/sf_jacobi.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_SF_JACOBI_TCC
#define _GLIBCXX_SF_JACOBI_TCC 1

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    _Tp
    __poly_jacobi(unsigned int __n, _Tp __alpha, _Tp __beta, _Tp __x)
    {
      auto _Pm2 = _Tp{1};
      if (__n == 0)
	return _Pm2;

      auto __apb = __alpha + __beta;
      auto __amb = __alpha - __beta;
      auto _Pm1 = (__amb + (_Tp{2} + __apb) * __x) / _Tp{2};
      if (__n == 1)
	return _Pm1;

      _Tp _P(0);
      for (auto __j = 2; __j < __n; ++__j )
	{
	  auto __japb = _Tp(__j) + __apb;
	  auto __japb1 = __japb + _Tp{1};
	  auto __japb2 = __japb1 + _Tp{1};
	  auto __c = _Tp{2} * _Tp(__j + 1) * __japb * __japb1;
	  auto __d = __japb1 * __apb * __amb;
	  auto __e = __japb * __japb1 * __japb2;
	  auto __f = _Tp{2} * (__j + __alpha) * (__j + __beta) * __japb2;

	  if (__c == _Tp{0})
	    std::__throw_runtime_error("poly_jacobi: error in recursion");
	  _P = ((__d + __e * __x) * _Pm1 - __f * _Pm2) / __c;
	  _Pm2 = _Pm1;
	  _Pm1 = _P;
	}
      return _P;
    }


  template<typename _Tp>
    _Tp
    __poly_radial_jacobi(unsigned int __n, unsigned int __m, _Tp __rho)
    {
      if (__m > __n)
	throw std::range_error("poly_radial_jacobi: m > n");
      else if ((__n - __m) % 2 == 1)
	return _Tp{0};
      else
	{
	  auto __k = (__n - __m) / 2;
	  return (__k % 2 == 0 ? 1 : -1) * std::pow(__rho, __m)
	       * __poly_jacobi(__k, _Tp(__m), _Tp{0},
			       _Tp{1} - _Tp{2} * __rho * __rho);
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_SF_JACOBI_TCC
