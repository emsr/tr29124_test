// Special functions -*- C++ -*-

// Copyright (C) 2006-2018 Free Software Foundation, Inc.
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

/** @file bits/specfun_util.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_SPECFUN_UTIL_H
#define _GLIBCXX_BITS_SPECFUN_UTIL_H 1

#pragma GCC system_header

#if __cplusplus >= 201103L
#  include <ratio>
#endif
#if __cplusplus >= 201402L
#  include <utility> // For exchange
#endif
#include <limits>

#define _GLIBCXX_HAVE_FLOAT128_MATH 0
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
#  if __has_include(<quadmath.h>)
#    include <quadmath.h>
#    define _GLIBCXX_HAVE_FLOAT128_MATH 1
#    if _GLIBCXX_USE_C99_MATH && !_GLIBCXX_USE_C99_FP_MACROS_DYNAMIC
namespace std
{
  bool isnan(__float128);
  bool isinf(__float128);
}
#    endif
#  endif
#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
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

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std


namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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

#if _GLIBCXX_HAVE_FLOAT128_MATH
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

// Need an exchange utility.
#if __cplusplus < 201402L
  template<typename _Tp, typename _Up = _Tp>
    _Tp
    __exchange(_Tp& __obj, _Up&& __new_val)
    {
      _Tp __old_val = std::move(__obj);
      __obj = std::forward<_Up>(__new_val);
      return __old_val;
    }
#else
  template<typename _Tp, typename _Up = _Tp>
    _Tp
    __exchange(_Tp& __obj, _Up&& __new_val)
    { return std::exchange(__obj, std::forward<_Up>(__new_val)); }
#endif

#endif // __cplusplus >= 201103L

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_SPECFUN_UTIL_H
