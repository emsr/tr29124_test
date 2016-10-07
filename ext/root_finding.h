// Math extensions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file ext/root_finding.h
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_ROOT_FINDING_H
#define _EXT_ROOT_FINDING_H 1

#pragma GCC system_header

#include <vector>
//#include <utility> // For pair?
#include <limits>

// Default __eps = _Tp{5} * std::numeric_limits<_Tp>::epsilon()?

// template<typename _Tp>
//   struct __bracket_t {_Tp __value_lower, _Tp __value_upper}; ?

// All the root methods should have the same API.
// brent, newton, safe don't have max_iter.

// __root_brackets should have an output iterator API?

// in __root_newton what if __df is zero?

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  // This struct contains bothe the value and the derivative at a point.
  template<typename _Tp>
    struct __root_state
    {
      _Tp __value;
      _Tp __deriv;
    };

  template<typename _Tp, typename _Func>
    bool
    __root_bracket(_Func __func, _Tp& __x_lower, _Tp& __x_upper,
    		   std::size_t __max_iter = 50);

  template<typename _Tp, typename _Func>
    std::vector<std::pair<_Tp, _Tp>>
    __root_brackets(_Func __func,
		    _Tp __x_lower, _Tp __x_upper, std::size_t __n);

  template<typename _Tp, typename _Func>
    _Tp
    __root_bisect(_Func __func, _Tp __x_lower, _Tp __x_upper, _Tp __eps,
		  std::size_t __max_iter = std::numeric_limits<_Tp>::digits);

  template<typename _Tp, typename _Func>
    _Tp
    __root_secant(_Func __func, _Tp __x_lower, _Tp __x_upper, _Tp __eps,
		  std::size_t __max_iter = 40);

  template<typename _Tp, typename _Func>
    _Tp
    __root_false_position(_Func __func, _Tp __x_lower, _Tp __x_upper,
			  _Tp __eps, std::size_t __max_iter = 40);

  template<typename _Tp, typename _Func>
    _Tp
    __root_ridder(_Func __func, _Tp __x_lower, _Tp __x_upper,
		  _Tp __eps, std::size_t __max_iter = 100);

  template<typename _Tp, typename _Func>
    _Tp
    __root_brent(_Func __func, _Tp __x_lower, _Tp __x_upper,
		 _Tp __eps, std::size_t __max_iter = 100);

  template<typename _Tp, typename _StateFunc>
    _Tp
    __root_newton(_StateFunc __func, _Tp __x_lower, _Tp __x_upper,
    		  _Tp __eps, std::size_t __max_iter = 40);

  template<typename _Tp, typename _StateFunc>
    _Tp
    __root_safe(_StateFunc __func, _Tp __x_lower, _Tp __x_upper,
    		_Tp __eps, std::size_t __max_iter = 100);

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <ext/root_finding.tcc>

#endif // _EXT_ROOT_FINDING_H
