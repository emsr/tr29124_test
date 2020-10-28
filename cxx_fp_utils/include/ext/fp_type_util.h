// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

/** @file ext/fp_type_util.h
 */

#ifndef _EXT_FP_TYPE_UTIL_H
#define _EXT_FP_TYPE_UTIL_H 1

#include <utility> // For exchange
#include <ratio>
#include <limits>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type.  This will be specialized for complex
   * and other types as appropriate.
   */
  template<typename _Tp>
    struct __num_traits
    {
      using __value_type = _Tp;
    };

  template<typename _Tp>
    using __num_traits_t = typename __num_traits<_Tp>::__value_type;


#if __cplusplus >= 201103L
  /**
   * Return a fraction as a real number.
   */
  template<intmax_t _Num, intmax_t _Den = 1, typename _Tp = double>
    inline constexpr _Tp
    __frac()
    {
      using __rat_t = std::ratio<_Num, _Den>;
      return _Tp(__rat_t::num) / _Tp(__rat_t::den);
    }
#endif


  /**
   * Create a NaN.
   * This will be overloaded for complex and vector types.
   */
  template<typename _Tp>
    struct __make_NaN
    {
      constexpr _Tp
      operator()()
      { return std::numeric_limits<_Tp>::quiet_NaN(); }
    };

#if __cplusplus >= 201103L

  /**
   * This is a more modern version of __promote_N in ext/type_traits.
   * This is used for numeric argument promotion of complex and cmath.
   */
  template<typename _Tp, bool = std::is_integral<_Tp>::value>
    struct fp_promote_help
    { using __type = double; };

  // No nested __type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename _Tp>
    struct fp_promote_help<_Tp, false>
    { };

  template<>
    struct fp_promote_help<float>
    { using __type = float; };

  template<>
    struct fp_promote_help<double>
    { using __type = double; };

  template<>
    struct fp_promote_help<long double>
    { using __type = long double; };

#ifdef _GLIBCXX_USE_FLOAT128
  template<>
    struct fp_promote_help<__float128>
    { using __type = __float128; };
#endif

  template<typename... _Tps>
    using fp_promote_help_t = typename fp_promote_help<_Tps...>::__type;

#if __cplusplus < 201402L
  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  template<typename _Tp, typename... _Tps>
    struct fp_promote
    {
      using __decaytype = typename std::decay<_Tp>::type;
      using __type = decltype(fp_promote_help_t<__decaytype>{}
		   + typename fp_promote<_Tps...>::__type{});
    };

  template<>
    template<typename _Tp>
      struct fp_promote<_Tp>
      {
	using __decaytype = typename std::decay<_Tp>::type;
	using __type = decltype(fp_promote_help_t<__decaytype>{});
      };
#else
  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  // @todo We use decay_t instead of remove_reference_t just in case _Tp
  // is an array.
  template<typename _Tp, typename... _Tps>
    struct fp_promote
    { using __type = decltype(fp_promote_help_t<std::decay_t<_Tp>>{}
		   + typename fp_promote<_Tps...>::__type{}); };

  template<typename _Tp>
    struct fp_promote<_Tp>
    { using __type = decltype(fp_promote_help_t<std::decay_t<_Tp>>{}); };
#endif

  template<typename... _Tps>
    using fp_promote_t = typename fp_promote<_Tps...>::__type;

#endif // __cplusplus >= 201103L

_GLIBCXX_END_NAMESPACE_VERSION
}

#endif // _EXT_FP_TYPE_UTIL_H
