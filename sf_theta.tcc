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

/** @file bits/sf_theta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_SF_THETA_TCC
#define _GLIBCXX_SF_THETA_TCC 1

#include <array>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    void
    __jacobi_cnsndn(_Tp __k, _Tp __u, _Tp& __cn, _Tp& __sn, _Tp& __dn)
    {
      _Tp __b, __c, __d;

      constexpr auto _S_N = 20;
      constexpr auto _S_CA = _Tp{0.00001};

      std::array<_Tp, _S_N> __m;
      std::array<_Tp, _S_N> __n;

      if (__k != _Tp{0})
	{
	  bool __bo = (__k < _Tp{0});
	  if (__bo)
	    {
	      __d = _Tp{1} - __k;
	      __k /= -_Tp{1} / __d;
	      __u *= (__d = std::sqrt(__d));
	    }
	  auto __a = _Tp{1};
	  __dn = _Tp{1};
	  auto __l = _S_N;
	  for (auto __i = 0; __i < _S_N; ++__i)
	    {
	      __l = __i;
	      __m[__i] = __a;
	      __n[__i] = (__k = std::sqrt(__k));
	      __c = 0.5 * (__a + __k);
	      if (std::abs(__a - __k) <= _S_CA * __a)
		break;
	      __k *= __a;
	      __a = __c;
	    }
	  __u *= __c;
	  __cn = std::cos(__u);
	  __sn = std::sin(__u);
	  if (__sn != _Tp{0})
	    {
	      __a = __cn / __sn;
	      __c *= __a;
	      for (auto __ii = __l; __ii + 1 >= 1; --__ii)
		{
		  _Tp __b = __m[__ii];
		  __a *= __c;
		  __c *= (__dn);
		  __dn = (__n[__ii] + __a) / (__b + __a);
		  __a = __c / __b;
		}
	      __a = _Tp{1} / std::sqrt(__c * __c + _Tp{1});
	      __sn = std::copysign(__a, __sn);
	      __cn = __c * __sn;
	    }
	  if (__bo)
	    {
	      __a = __dn;
	      __dn = __cn;
	      __cn = __a;
	      __sn /= __d;
	    }
	}
      else
	{
	  __cn = _Tp{1} / std::cosh(__u);
	  __sn = std::tanh(__u);
	  __dn = __cn;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_SF_THETA_TCC
