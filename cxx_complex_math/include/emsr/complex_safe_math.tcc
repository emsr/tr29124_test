// Special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

/** @file ext/complex_safe_math.tcc
 */

#ifndef COMPLEX_SAFE_MATH_TCC
#define COMPLEX_SAFE_MATH_TCC 1

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{

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
  template<typename Tp>
    std::complex<Tp>
    safe_div(const std::complex<Tp>& z1, const std::complex<Tp>& z2)
    {
      // Half the largest available floating-point number.
      const auto s_hmax = emsr::lim_max(std::real(z1)) / Tp{2};

      auto re1 = std::real(z1);
      auto im1 = std::imag(z1);
      auto re2 = std::real(z2);
      auto im2 = std::imag(z2);

      //  Find the largest and smallest magnitudes
      auto z1b = std::max(std::abs(re1), std::abs(im1));
      auto z2max = std::abs(re2);
      auto z2min = std::abs(im2);
      if (z2max < z2min)
	std::swap(z2max, z2min);

      if (z2max < Tp{1} && z1b > z2max * s_hmax)
	throw std::runtime_error(__N("safe_div: overflow in complex division"));

      re1 /= z1b;
      im1 /= z1b;
      re2 /= z2max;
      im2 /= z2max;
      auto term = z2min / z2max;
      auto denom = Tp{1} + term * term;
      auto scale = z1b / z2max / denom;
      auto qr = (re1 * re2 + im1 * im2) * scale;
      auto qi = (re2 * im1 - re1 * im2) * scale;

      return std::complex<Tp>{qr, qi};
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
  template<typename Tp>
    Tp
    safe_mul(Tp s1, Tp s2)
    {
      // The largest available floating-point number.
      const auto s_max = emsr::lim_max(std::real(s1));
      const auto s_sqrt_max = emsr::sqrt_max(std::real(s1));
      auto abs_s1 = std::abs(s1);
      auto abs_s2 = std::abs(s2);
      if (abs_s1 < s_sqrt_max || abs_s2 < s_sqrt_max)
	{
	  auto abs_max = abs_s1;
	  auto abs_min = abs_s2;
	  if (abs_max < abs_min)
	    std::swap(abs_max, abs_min);
	  if (abs_max > s_sqrt_max && abs_min > s_max / abs_max)
	    throw std::runtime_error(__N("safe_mul: overflow in scalar multiplication"));
	  else
	    return s1 * s2;
	}
      else
	throw std::runtime_error(__N("safe_mul: overflow in scalar multiplication"));
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
  template<typename Tp>
    std::complex<Tp>
    safe_mul(const std::complex<Tp>& z1, const std::complex<Tp>& z2)
    {
      // Half the largest available floating-point number.
      const auto s_max = emsr::lim_max(std::real(z1));
      const auto s_sqrt_max = emsr::sqrt_max(std::real(z1));

      auto re1 = std::real(z1);
      auto im1 = std::imag(z1);
      auto re2 = std::real(z2);
      auto im2 = std::imag(z2);

      auto abs_rem = std::abs(re1 - im1);
      auto abs_rep = std::abs(re2 + im2);
      if (abs_rem < s_sqrt_max || abs_rep < s_sqrt_max)
	{
	  // Find the largest and smallest magnitudes
	  auto abs_min = abs_rem;
	  auto abs_max = abs_rep;
	  if (abs_max < abs_min)
	    std::swap(abs_max, abs_min);
	  if (abs_max > s_sqrt_max && abs_min > s_max / abs_max)
	    throw std::runtime_error(__N("safe_mul: overflow in complex multiplication"));
	  else
	    return std::complex<Tp>((re1 - im1) * (re2 + im2),
			safe_mul(re1, im2) + safe_mul(re2, im1));
	}
      else
	throw std::runtime_error(__N("safe_mul: overflow in complex multiplication"));
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
  template<typename Tp>
    std::complex<Tp>
    safe_sqr(const std::complex<Tp>& z)
    {
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_max = emsr::lim_max<Tp>();
      const auto s_hmax = s_max / Tp{2};
      const auto s_sqrt_max = emsr::sqrt_max<Tp>();
      const auto s_sqrt_hmax = s_sqrt_max / s_sqrt_2;

      auto rez = std::real(z);
      auto imz = std::imag(z);
      auto abs_rez = std::abs(rez);
      auto abs_imz = std::abs(imz);
      auto zm = rez - imz;
      auto zp = rez + imz;
      auto abs_zm = std::abs(zm);
      auto abs_zp = std::abs(zp);

      if ((abs_zm < s_sqrt_max || abs_zp < s_sqrt_max)
       && (abs_rez < s_sqrt_hmax || abs_imz < s_sqrt_hmax))
	{
	  // Sort the magnitudes of the imag part factors.
	  auto imzmax = abs_rez;
	  auto imzmin = abs_imz;
	  if (imzmax < imzmin)
	    std::swap(imzmax, imzmin);
	  if (imzmax >= s_sqrt_hmax && imzmin > s_hmax / imzmax)
	    throw std::runtime_error(__N("safe_sqr: overflow in complex multiplication"));

	  // Sort the magnitudes of the real part factors.
	  auto rezmax = abs_zp;
	  auto rezmin = abs_zm;
	  if (imzmax < rezmin)
	    std::swap(rezmax, rezmin);
	  if (rezmax >= s_sqrt_max && rezmin > s_max / rezmax)
	    throw std::runtime_error(__N("safe_sqr: "
					"overflow in complex multiplication"));

	  return std::complex<Tp>(zm * zp, Tp{2} * rez * imz);
	}
      else
	throw std::runtime_error(__N("safe_sqr: overflow in complex multiplication"));
    }

} // namespace emsr

#endif // COMPLEX_SAFE_MATH_TCC
