
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

#ifndef FP_TYPE_UTIL_H
#define FP_TYPE_UTIL_H 1

#include <utility> // For exchange
#include <ratio>
#include <limits>

namespace emsr
{

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type.  This will be specialized for complex
   * and other types as appropriate.
   */
  template<typename Tp>
    struct num_traits
    {
      using value_type = Tp;
    };

  template<typename Tp>
    using num_traits_t = typename num_traits<Tp>::value_type;


  /**
   * Return a fraction as a real number.
   */
  template<intmax_t Num, intmax_t Den = 1, typename Tp = double>
    inline constexpr Tp
    frac()
    {
      using rat_t = std::ratio<Num, Den>;
      return Tp(rat_t::num) / Tp(rat_t::den);
    }


  /**
   * Create a NaN.
   * This will be overloaded for complex and vector types.
   */
  template<typename Tp>
    struct make_NaN
    {
      constexpr Tp
      operator()()
      { return std::numeric_limits<Tp>::quiet_NaN(); }
    };

  /**
   * This is a more modern version of promote_N in ext/type_traits.
   * This is used for numeric argument promotion of complex and cmath.
   */
  template<typename Tp, bool = std::is_integral_v<Tp>>
    struct fp_promote_help
    { using type = double; };

  // No nested type member for non-integer non-floating point types,
  // allows this type to be used for SFINAE to constrain overloads in
  // <cmath> and <complex> to only the intended types.
  template<typename Tp>
    struct fp_promote_help<Tp, false>
    { }; // using type = void?

  template<>
    struct fp_promote_help<float>
    { using type = float; };

  template<>
    struct fp_promote_help<double>
    { using type = double; };

  template<>
    struct fp_promote_help<long double>
    { using type = long double; };

#ifdef EMSR_HAVE_FLOAT128
  template<>
    struct fp_promote_help<__float128>
    { using type = __float128; };
#endif

  template<typename... _Tps>
    using fp_promote_help_t = typename fp_promote_help<_Tps...>::type;

  // Decay refs and cv...
  // Alternatively we could decay refs and propagate cv to promoted type.
  // @todo We use decay_t instead of remove_reference_t just in case Tp
  // is an array.
  template<typename Tp, typename... _Tps>
    struct fp_promote
    {
      using type = decltype(fp_promote_help_t<std::decay_t<Tp>>{}
		          + typename fp_promote<_Tps...>::type{});
    };

  template<typename Tp>
    struct fp_promote<Tp>
    {
      using type = decltype(fp_promote_help_t<std::decay_t<Tp>>{});
    };

  template<typename... _Tps>
    using fp_promote_t = typename fp_promote<_Tps...>::type;

}

#endif // FP_TYPE_UTIL_H
