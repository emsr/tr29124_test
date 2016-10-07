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

/** @file bits/complex_util.tcc
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_COMPLEX_UTIL_TCC
#define _GLIBCXX_BITS_COMPLEX_UTIL_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief Carefully compute and return @c z1/z2 avoiding overflow
   *       and destructive underflow.
   *	   If the quotient can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   *
   * @param[in]  z1  Dividend
   * @param[in]  z2  Divisor
   * @return  The quotient of z1 and z2
   * @throws  std::runtime_error on division overflow.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_div(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2)
    {
      // Half the largest available floating-point number.
      constexpr _Tp _S_hmax = __gnu_cxx::__max<_Tp>() / _Tp{2};

      auto __re1 = std::real(__z1);
      auto __im1 = std::imag(__z1);
      auto __re2 = std::real(__z2);
      auto __im2 = std::imag(__z2);

      //  Find the largest and smallest magnitudes
      auto __z1b = std::max(std::abs(__re1), std::abs(__im1));
      auto __z2max = std::abs(__re2);
      auto __z2min = std::abs(__im2);
      if (__z2max < __z2min)
	std::swap(__z2max, __z2min);

      if (__z2max < _Tp{1} && __z1b > __z2max * _S_hmax)
	std::__throw_runtime_error(__N("__safe_div: "
				       "overflow in complex division"));

      __re1 /= __z1b;
      __im1 /= __z1b;
      __re2 /= __z2max;
      __im2 /= __z2max;
      auto __term = __z2min / __z2max;
      auto __denom = _Tp{1} + __term * __term;
      auto __scale = __z1b / __z2max / __denom;
      auto __qr = (__re1 * __re2 + __im1 * __im2) * __scale;
      auto __qi = (__re2 * __im1 - __re1 * __im2) * __scale;

      return std::complex<_Tp>{__qr, __qi};
    }

  /**
   * @brief Carefully compute and return @c s1*s2 avoiding overflow.
   *	   If the product can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   *
   * @param[in]  s1  Factor 1
   * @param[in]  s2  Factor 2
   * @return  The product of s1 and s2
   * @throws  std::runtime_error on multiplication overflow.
   */
  template<typename _Tp>
    _Tp
    __safe_mul(_Tp __s1, _Tp __s2)
    {
      // The largest available floating-point number.
      const _Tp _S_max = __gnu_cxx::__max<_Tp>();
      const auto _S_sqrt_max = __gnu_cxx::__sqrt_max<_Tp>();
      auto __abs_s1 = std::abs(__s1);
      auto __abs_s2 = std::abs(__s2);
      if (__abs_s1 < _S_sqrt_max || __abs_s2 < _S_sqrt_max)
	{
	  auto __abs_max = __abs_s1;
	  auto __abs_min = __abs_s2;
	  if (__abs_max < __abs_min)
	    std::swap(__abs_max, __abs_min);
	  if (__abs_max > _S_sqrt_max && __abs_min > _S_max / __abs_max)
	    std::__throw_runtime_error(__N("__safe_mul: "
					"overflow in scalar multiplication"));
	  else
	    return __s1 * __s2;
	}
      else
	std::__throw_runtime_error(__N("__safe_mul: "
				       "overflow in scalar multiplication"));
    }

  /**
   * @brief Carefully compute and return @c z1*z2 avoiding overflow.
   *	   If the product can be successfully computed it is returned.
   *	   Otherwise, std::runtime_error is thrown.
   *
   * @param[in]  z1  Factor 1
   * @param[in]  z2  Factor 2
   * @return  The product of z1 and z2
   * @throws  std::runtime_error on multiplication overflow.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_mul(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2)
    {
      // Half the largest available floating-point number.
      const _Tp _S_max = __gnu_cxx::__max<_Tp>();
      const auto _S_sqrt_max = __gnu_cxx::__sqrt_max<_Tp>();

      auto __re1 = std::real(__z1);
      auto __im1 = std::imag(__z1);
      auto __re2 = std::real(__z2);
      auto __im2 = std::imag(__z2);

      auto __abs_rem = std::abs(__re1 - __im1);
      auto __abs_rep = std::abs(__re2 + __im2);
      if (__abs_rem < _S_sqrt_max || __abs_rep < _S_sqrt_max)
	{
	  // Find the largest and smallest magnitudes
	  auto __abs_min = __abs_rem;
	  auto __abs_max = __abs_rep;
	  if (__abs_max < __abs_min)
	    std::swap(__abs_max, __abs_min);
	  if (__abs_max > _S_sqrt_max && __abs_min > _S_max / __abs_max)
	    std::__throw_runtime_error(__N("__safe_mul: "
					"overflow in complex multiplication"));
	  else
	    return std::complex<_Tp>((__re1 - __im1) * (__re2 + __im2),
			__safe_mul(__re1, __im2) + __safe_mul(__re2, __im1));
	}
      else
	std::__throw_runtime_error(__N("__safe_mul: "
				       "overflow in complex multiplication"));
    }

  /**
   * @brief Carefully compute @c z*z avoiding overflow.
   *	    If the product can be successfully computed it is returned.
   *	    Otherwise, std::runtime_error is thrown.
   *
   * @param[in]  z  Argument
   * @return  The square of the argument
   * @throws  std::runtime_error on multiplication overflow.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __safe_sqr(const std::complex<_Tp>& __z)
    {
      const auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      const auto _S_max = __gnu_cxx::__max<_Tp>();
      const auto _S_hmax = _S_max / _Tp{2};
      const auto _S_sqrt_max = __gnu_cxx::__sqrt_max<_Tp>();
      const auto _S_sqrt_hmax = _S_sqrt_max / _S_sqrt_2;

      auto __rez = std::real(__z);
      auto __imz = std::imag(__z);
      auto __abs_rez = std::abs(__rez);
      auto __abs_imz = std::abs(__imz);
      auto __zm = __rez - __imz;
      auto __zp = __rez + __imz;
      auto __abs_zm = std::abs(__zm);
      auto __abs_zp = std::abs(__zp);

      if ((__abs_zm < _S_sqrt_max || __abs_zp < _S_sqrt_max)
       && (__abs_rez < _S_sqrt_hmax || __abs_imz < _S_sqrt_hmax))
	{
	  // Sort the magnitudes of the imag part factors.
	  auto __imzmax = __abs_rez;
	  auto __imzmin = __abs_imz;
	  if (__imzmax < __imzmin)
	    std::swap(__imzmax, __imzmin);
	  if (__imzmax >= _S_sqrt_hmax && __imzmin > _S_hmax / __imzmax)
	    std::__throw_runtime_error(__N("__safe_sqr: "
					"overflow in complex multiplication"));

	  // Sort the magnitudes of the real part factors.
	  auto __rezmax = __abs_zp;
	  auto __rezmin = __abs_zm;
	  if (__imzmax < __rezmin)
	    std::swap(__rezmax, __rezmin);
	  if (__rezmax >= _S_sqrt_max && __rezmin > _S_max / __rezmax)
	    std::__throw_runtime_error(__N("__safe_sqr: "
					"overflow in complex multiplication"));

	  return std::complex<_Tp>(__zm * __zp, _Tp{2} * __rez * __imz);
	}
      else
	std::__throw_runtime_error(__N("__safe_sqr: "
				       "overflow in complex multiplication"));
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_COMPLEX_UTIL_TCC
