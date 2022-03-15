
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

/** @file bits/sf_airy.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef SF_CARDINAL_TCC
#define SF_CARDINAL_TCC 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_trig.h>

namespace emsr
{
namespace detail
{

  /**
   * @brief Return the sinus cardinal function
   * @f[
   *   sinc(x) = \frac{\sin(x)}{x}
   * @f]
   */
  template<typename Tp>
    emsr::fp_promote_t<Tp>
    sinc(Tp x)
    {
      if (std::isnan(x))
        return emsr::quiet_NaN(x);
      else if (std::abs(x) == emsr::infinity(x))
	return Tp{0};
      else if (std::abs(x) < emsr::sqrt_min(x))
        return Tp{1} - x * x / Tp{6};
      else
	return std::sin(x) / x;
    }

  /**
   * @brief Return the reperiodized sinus cardinal function
   * @f[
   *   sinc_\pi(x) = \frac{\sin(\pi x)}{\pi x}
   * @f]
   */
  template<typename Tp>
    emsr::fp_promote_t<Tp>
    sinc_pi(Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      if (std::isnan(x))
        return emsr::quiet_NaN(x);
      else if (std::abs(x) == emsr::infinity(x))
	return Tp{0};
      else
	{
	  auto arg = s_pi * x;
	  if (std::abs(arg) < Tp{4} * emsr::sqrt_min(x))
	    return Tp{1} - arg * arg / Tp{6};
	  else
	    return emsr::sin_pi(x) / arg;
	}
    }

  /**
   * @brief Return the hyperbolic sinus cardinal function
   * @f[
   *   sinhc(x) = \frac{\sinh(x)}{x}
   * @f]
   */
  template<typename Tp>
    emsr::fp_promote_t<Tp>
    sinhc(Tp x)
    {
      if (std::isnan(x))
        return emsr::quiet_NaN(x);
      else if (std::abs(x) < Tp{4} * emsr::sqrt_min(x))
        return Tp{1} + x * x / Tp{6};
      else
	return std::sinh(x) / x;
    }

  /**
   * @brief Return the reperiodized hyperbolic sinus cardinal function
   * @f[
   *   sinhc_\pi(x) = \frac{\sinh(\pi x)}{\pi x}
   * @f]
   */
  template<typename Tp>
    emsr::fp_promote_t<Tp>
    sinhc_pi(Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      if (std::isnan(x))
        return emsr::quiet_NaN(x);
      else
	{
	  auto arg = s_pi * x;
	  if (std::abs(arg) < Tp{4} * emsr::sqrt_min(x))
	    return Tp{1} + arg * arg / Tp{6};
	  else
	    return emsr::sinh_pi(x) / arg;
	}
    }
} // namespace detail
} // namespace emsr

#endif // SF_CARDINAL_TCC
