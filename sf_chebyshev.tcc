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

/** @file bits/sf_chebyshev.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_SF_CHEBYSHEV_TCC
#define _GLIBCXX_SF_CHEBYSHEV_TCC 1

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    _Tp
    __chebyshev_recur(unsigned int __n, _Tp __x, _Tp _C0, _Tp _C1)
    {
      auto _C = _Tp(0);
      for (unsigned int __j = 1; __j < __n; ++__j)
      {
	_C = _Tp(2) * __x * _C1 - _C0;
	_C0 = _C1;
	_C1 = _C;
      }
      return _C;
    }

  template<typename _Tp>
    _Tp
    __chebyshev_t(unsigned int __n, _Tp __x)
    {
      auto _T0 = _Tp(1);
      if (__n == 0)
	return _T0;

      auto _T1 = __x;
      if (__n == 1)
	return _T1;

      return __chebyshev_recur(__n, __x, _T0, _T1);
    }

  template<typename _Tp>
    _Tp
    __chebyshev_u(unsigned int __n, _Tp __x)
    {
      auto _U0 = _Tp(1);
      if (__n == 0)
	return _U0;

      auto _U1 = _Tp(2) * __x;
      if (__n == 1)
	return _U1;

      return __chebyshev_recur(__n, __x, _U0, _U1);
    }

  template<typename _Tp>
    _Tp
    __chebyshev_v(unsigned int __n, _Tp __x)
    {
      auto _V0 = _Tp(1);
      if (__n == 0)
	return _V0;

      auto _V1 = _Tp(2) * __x - _Tp(1);
      if (__n == 1)
	return _V1;

      return __chebyshev_recur(__n, __x, _V0, _V1);
    }

  template<typename _Tp>
    _Tp
    __chebyshev_w(unsigned int __n, _Tp __x)
    {
      auto _W0 = _Tp(1);
      if (__n == 0)
	return _W0;

      auto _W1 = _Tp(2) * __x + _Tp(1);
      if (__n == 1)
	return _W1;

      return __chebyshev_recur(__n, __x, _W0, _W1);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
}

#endif // _GLIBCXX_SF_CHEBYSHEV_TCC
