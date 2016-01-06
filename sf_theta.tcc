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

/** @file bits/sf_theta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_SF_THETA_TCC
#define _GLIBCXX_SF_THETA_TCC 1

#include <vector>
#include <tuple>
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    std::tuple<_Tp, _Tp, _Tp>
    __jacobi_sncndn(_Tp __k, _Tp __u)
    {
      using _Val = __num_traits_t<_Tp>;
      constexpr auto _S_eps = __gnu_cxx::__math_constants<_Val>::__eps;
      constexpr auto _S_NaN = std::numeric_limits<_Val>::quiet_NaN();

      if (__isnan(__k) || __isnan(__u))
	return std::make_tuple(_S_NaN, _S_NaN, _S_NaN);
      else if (std::abs(__k) > _Tp{1})
	throw std::domain_error("__jacobi_sncndn:"
			" argument k out of range");
      else if (std::abs(_Tp{1} - __k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::tanh(__u);
	  auto __cn = _Tp{1} / std::cosh(__u);
	  auto __dn = __cn;
	  return std::make_tuple(__sn, __cn, __dn);
	}
      else if (std::abs(__k) < _Tp{2} * _S_eps)
	{
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
	  auto __dn = _Tp{1};
	  return std::make_tuple(__sn, __cn, __dn);
	}
      else
	{
	  constexpr auto _S_CA = std::sqrt(_S_eps);
	  constexpr auto _S_N = 100;
	  std::vector<_Tp> __m;
	  std::vector<_Tp> __n;
	  __m.reserve(20);
	  __n.reserve(20);
	  _Tp __c, __d;
	  auto __mc = _Tp{1} - __k * __k;
	  bool __bo = (__mc < _Tp{0});
	  if (__bo)
	    {
	      __d = _Tp{1} - __mc;
	      __mc /= -_Tp{1} / __d;
	      __u *= (__d = std::sqrt(__d));
	    }
	  auto __a = _Tp{1};
	  auto __dn = _Tp{1};
	  auto __l = _S_N;
	  for (auto __i = 0; __i < _S_N; ++__i)
	    {
	      __l = __i;
	      __m.push_back(__a);
	      __n.push_back(__mc = std::sqrt(__mc));
	      __c = 0.5 * (__a + __mc);
	      if (std::abs(__a - __mc) <= _S_CA * __a)
		break;
	      __mc *= __a;
	      __a = __c;
	    }
	  __u *= __c;
	  auto __sn = std::sin(__u);
	  auto __cn = std::cos(__u);
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
	      __a = _Tp{1} / std::hypot(_Tp{1}, __c);
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
	  return std::make_tuple(__sn, __cn, __dn);
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_SF_THETA_TCC
