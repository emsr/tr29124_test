// TR29124 math special functions -*- C++ -*-

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

/** @file bits/complex_util.h
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_COMPLEX_UTIL_H
#define _GLIBCXX_BITS_COMPLEX_UTIL_H 1

#pragma GCC system_header

#include <complex>
#include <ratio>
#include <limits>
#include <bits/specfun_util.h>
#include <bits/numeric_limits.h>
#include <ext/math_util.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

//
// The base definitions of __num_traits, __promote_fp,
// and __isnan reside in ext/specfun_util.h
//

  /**
   * A class to reach into compound numeric types to extract the
   * value or element type - specialized for complex.
   */
  template<>
    template<typename _Tp>
      struct __num_traits<std::complex<_Tp>>
      {
	using __value_type = typename std::complex<_Tp>::value_type;
      };

  /**
   * Create a complex number NaN.
   */
  template<>
    template<typename _Tp>
      struct __make_NaN<std::complex<_Tp>>
      {
	constexpr std::complex<_Tp>
	operator()()
	{
	  auto __NaN = std::numeric_limits<_Tp>::quiet_NaN();
	  return std::complex<_Tp>{__NaN, __NaN};
	}
      };

  /**
   * Return true if one component of a complex number is NaN.
   */
  template<typename _Tp>
    inline bool
    __isnan(const std::complex<_Tp>& __z)
    { return __isnan(std::real(__z)) || __isnan(std::imag(__z)); }

  /**
   * Return true if one component of a complex number is inf.
   */
  template<typename _Tp>
    inline bool
    __isinf(const std::complex<_Tp>& __z)
    { return __isinf(std::real(__z)) || __isinf(std::imag(__z)); }


  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l1_norm(const std::complex<_Tp>& __z)
    { return std::abs(std::real(__z)) + std::abs(std::imag(__z)); }

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l2_norm(const std::complex<_Tp>& __z)
    { return std::norm(__z); }

  /**
   * Return the Linf norm modulus of a complex number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __linf_norm(const std::complex<_Tp>& __z)
    { return std::max(std::abs(std::real(__z)), std::abs(std::imag(__z))); }


  /**
   * Return the L1 norm modulus or the Manhattan metric distance of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l1_norm(_Tp __x)
    { return std::abs(__x); }

  /**
   * Return the L2 norm modulus or the Euclidean metric distance of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __l2_norm(_Tp __x)
    { return std::abs(__x); }

  /**
   * Return the Linf norm modulus of a real number.
   */
  template<typename _Tp>
    inline constexpr _Tp
    __linf_norm(_Tp __x)
    { return std::abs(__x); }


  /**
   * Carefully compute @c z1/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computedit is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_div(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2);

  /**
   * Carefully compute @c s/z2 avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_div(_Sp __s, const std::complex<_Tp>& __z)
    { return __safe_div(std::complex<_Tp>(__s), __z); }

  /**
   * Carefully compute @c z1/s avoiding overflow and destructive underflow.
   * If the quotient is successfully computed it is returned.
   * Otherwise, @c false is returned and the quotient is not.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_div(const std::complex<_Tp>& __z, _Sp __s)
    { return __safe_div(__z, std::complex<_Tp>(__s)); }

  /**
   * @brief Carefully compute and return @c s1*s2 avoiding overflow.
   *	   If the product can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    _Tp
    __safe_mul(_Tp __s1, _Tp __s2);

  /**
   * Carefully compute @c z1*z2 avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_mul(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2);

  /**
   * Carefully compute @c s*z avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_mul(_Sp __s, const std::complex<_Tp>& __z)
    { return __safe_mul(std::complex<_Tp>(__s), __z); }

  /**
   * Carefully compute @c z*s avoiding overflow.
   * If the product is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Sp, typename _Tp>
    inline std::complex<_Tp>
    __safe_mul(const std::complex<_Tp>& __z, _Sp __s)
    { return __safe_mul(__z, std::complex<_Tp>(__s)); }

  /**
   * Carefully compute @c z*z avoiding overflow.
   * If the square is successfully computed it is returned.
   * Otherwise, std::runtime_error is thrown.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_sqr(const std::complex<_Tp>& __z);


  /**
   * A function to reliably test a complex number for realness.
   *
   * @param __w  The complex argument
   * @param __mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Im(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename _Tp>
    bool
    __fpreal(const std::complex<_Tp>& __w, const _Tp __mul = _Tp{5})
    { return __gnu_cxx::__fpequal(std::imag(__w), _Tp{0}, __mul); }

  // Specialize for real numbers.
  template<typename _Tp>
    bool
    __fpreal(const _Tp)
    { return true; }


  /**
   * A function to reliably test a complex number for imaginaryness [?].
   *
   * @param __w  The complex argument
   * @param __mul The multiplier for numeric epsilon for comparison tolerance
   * @return  @c true if @f$ Re(w) @f$ is zero within @f$ mul * epsilon @f$,
   *          @c false otherwize.
   */
  template<typename _Tp>
    bool
    __fpimag(const std::complex<_Tp>& __w, const _Tp __mul = _Tp{5})
    { return __gnu_cxx::__fpequal(std::real(__w), _Tp{0}, __mul); }

  // Specialize for real numbers.
  template<typename _Tp>
    bool
    __fpimag(const _Tp)
    { return false; }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

#if __cplusplus >= 201103L

  /**
   * We need isnan to be extended to std::complex.
   */
  template<typename _Tp>
    bool
    isnan(const std::complex<_Tp>& __z)
    { return std::__detail::__isnan(__z); }

  /**
   * This is a more modern version of __promote_N in ext/type_traits
   * specialized for complex.
   * This is used for numeric argument promotion of complex and cmath
   */
  template<>
    template<typename _Tp>
      struct __promote_fp_help<std::complex<_Tp>, false>
      {
      private:
	using __vtype = typename std::complex<_Tp>::value_type;
      public:
	using __type = decltype(std::complex<__promote_fp_help_t<__vtype>>{});
      };

  /**
   * Type introspection for complex.
   */
  template<typename _Tp>
    struct is_complex : public std::false_type
    { };

  /**
   * Type introspection for complex.
   */
  template<>
    template<typename _Tp>
      struct is_complex<std::complex<_Tp>> : public std::true_type
      { };

  /**
   * Type introspection for complex.
   */
  template<typename _Tp>
    using is_complex_t = typename is_complex<_Tp>::type;

  /**
   * Type introspection for complex.
   */
  template<typename _Tp>
    constexpr bool is_complex_v = is_complex<_Tp>::value;

#endif // __cplusplus >= 201103L

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <bits/complex_util.tcc>

#endif // _GLIBCXX_BITS_COMPLEX_UTIL_H
