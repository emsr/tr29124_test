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

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief Carefully compute @c z1/z2 avoiding overflow
   *       and destructive underflow.
   *	   If the quotient can be successfully computed, @c true
   *	   is returned and the quotient is returned in @c z1dz2.
   *	   Otherwise, std::runtime_error is thrown.
   *
   * @todo  In a better world optional might be a batter solution.
   *        This would be nothrow...
   *
   * @param[in]  z1  Dividend
   * @param[in]  z2  Divisor
   * @param[out]  z1dz2  Quotient
   * @return  @c true on success
   * @throws  std::runtime_error on division overflow.
   */
  template<typename _Tp>
    bool
    __safe_div(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2,
	       std::complex<_Tp>& __z1dz2)
    {
      //  Note that _S_xhinf is a machine floating-point dependent constant
      //  set equal to half the largest available floating-point number.
      static constexpr _Tp _S_xhinf = _Tp(0.5) * std::numeric_limits<_Tp>::max();

      auto __re1 = std::real(__z1);
      auto __im1 = std::imag(__z1);
      auto __re2 = std::real(__z2);
      auto __im2 = std::imag(__z2);

      //  Find the largest and smallest magnitudes
      auto __z1b = std::max(std::abs(__re1), std::abs(__im1));
      //std::tie<auto>(__z2min, __z2max) = std::minmax(std::abs(__re2), std::abs(__im2));
      auto __z2max = std::abs(__re2);
      auto __z2min = std::abs(__im2);

      if (__z2max < __z2min)
	std::swap(__z2max, __z2min);

      if (__z2max < _Tp(1) && __z1b > __z2max * _S_xhinf)
	std::__throw_runtime_error(__N("__safe_div: "
				       "overflow in complex division"));

      __re1 /= __z1b;
      __im1 /= __z1b;
      __re2 /= __z2max;
      __im2 /= __z2max;
      auto __term = __z2min / __z2max;
      auto __denom = _Tp(1) + __term * __term;
      auto __scale = __z1b / __z2max / __denom;
      auto __qr = (__re1 * __re2 + __im1 * __im2) * __scale;
      auto __qi = (__re2 * __im1 - __re1 * __im2) * __scale;
      __z1dz2 = std::complex<_Tp>{__qr, __qi};

      return true;
    }

  /**
   * @brief Carefully compute @c z1*z2 avoiding overflow
   *       and destructive underflow.
   *	   If the product can be successfully computed, @c true
   *	   is returned and the product is returned in @c z1z2.
   *	   Otherwise, std::runtime_error is thrown.
   *
   * @todo  In a better world optional might be a batter solution.
   *        This would be nothrow...
   *
   * @param[in]  z1  Mutiplicand 1
   * @param[in]  z2  Mutiplicand 2
   * @param[out]  z1z2  Product
   * @return  @c true on success
   * @throws  std::runtime_error on multiplication overflow.
   */
  template<typename _Tp>
    bool
    __safe_mul(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2,
	       std::complex<_Tp>& __z1z2)
    {
      //  Note that _S_xhinf is a machine floating-point dependent constant
      //  set equal to half the largest available floating-point number.
      static constexpr _Tp _S_xhinf = _Tp(0.5) * std::numeric_limits<_Tp>::max();

      auto __re1 = std::real(__z1);
      auto __im1 = std::imag(__z1);
      auto __re2 = std::real(__z2);
      auto __im2 = std::imag(__z2);

      if (_abs_re1 < _S_sqrt_max && _abs_Im_zhat < _S_sqrt_max)
	{
      //  Find the largest and smallest magnitudes
      auto __z1b = std::max(std::abs(__re1), std::abs(__im1));
      auto __z2max = std::abs(__re2);
      auto __z2min = std::abs(__im2);

      if (__z2max < __z2min)
	std::swap(__z2max, __z2min);


      return true;
    }

  /**
   * @brief Carefully compute @c z*z avoiding overflow and destructive underflow.
   *	    If the product can be successfully computed, @c true
   *	    is returned and the square is returned in @c z2.
   *	    Otherwise, std::runtime_error is thrown.
   *
   * @todo  In a better world optional might be a batter solution.
   *        This would be nothrow...
   *
   * @param[in]  z  Argument
   * @param[out]  z2  Square of the argument
   * @return  @c true on success
   * @throws  std::runtime_error on multiplication overflow.
   */
  template<typename _Tp>
    bool
    __safe_sqr(const std::complex<_Tp>& __z1, const std::complex<_Tp>& __z2,
	       std::complex<_Tp>& __z2)
    {
      //  Note that _S_xhinf is a machine floating-point dependent constant
      //  set equal to half the largest available floating-point number.
      static constexpr _Tp _S_xhinf = _Tp(0.5) * std::numeric_limits<_Tp>::max();

      auto __rez = std::real(__z);
      auto __imz = std::imag(__z);
      auto __abs_rez = std::abs(__re);
      auto __abs_imz = std::abs(__im);

      if (_abs_rez < _S_sqrt_max && _abs_imz < _S_sqrt_max)
	{
	  //  Find the largest and smallest magnitudes
	  auto __zmax = __abs_rez;
	  auto __zmin = __abs_imz;
	  if (__zmax < __zmin)
	    std::swap(__zmax, __zmin);
	  //std::tie<auto>(__zmin, __zmax) = std::minmax(__abs_rez, __abs_imz);
	  if (__zmax >= _S_1d2 && __zmin > _S_inf / (_Tp{2} * __zmax))
	    std::__throw_runtime_error(__N("__safe_sqr: "
					 "overflow in complex multiplication"));
	  else
	    __z2 = std::complex<_Tp>((__rez - __imz) * (__rez + __imz),
				     _Tp{2} * __rez * __imz);
	}
      else
	std::__throw_runtime_error(__N("__safe_sqr: "
				       "overflow in complex multiplication"));

      return true;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_COMPLEX_UTIL_TCC
