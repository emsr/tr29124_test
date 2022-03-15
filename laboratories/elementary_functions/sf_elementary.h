
// Copyright (C) 2022 Edward M. Smith-Rowland
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

#ifndef SF_ELEMENTARY_H
#define SF_ELEMENTARY_H 1

#include <cmath>
#include <complex>

#include "sf_elementary.tcc"

namespace emsr
{

/**
 * Return the Wright omega function for complex argument.
 */
template<typename Tp>
  std::complex<Tp>
  wright_omega(const std::complex<Tp>& z)
  {
    std::complex<Tp> err, res, condn;
    return detail::wright_omega(z, err, res, condn);
  }

/**
 * Return the Wright omega function for real argument.
 */
template<typename Tp>
  Tp
  wright_omega(Tp x)
  {
    using Cmplx = std::complex<Tp>;
    Cmplx z(x), err, res, condn;
    return std::real(detail::wright_omega(z, err, res, condn));
  }

/**
 * Return the 0-branch of the Lambert W function of real argument.
 * This is denoted Wp in the DLMF.
 */
template<typename Tp>
  Tp
  lambert_wp(Tp x)
  {
    using Cmplx = std::complex<Tp>;
    const auto s_i = Cmplx(Tp{0}, Tp{1});
    const auto s_1de = emsr::inv_e_v<Tp>;
    if (x < -s_1de)
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (x == -s_1de)
      return Tp{-1};
    else
      {
	Cmplx z = std::log(x + s_i * std::numeric_limits<Tp>::epsilon());
	Cmplx err, res, condn;
	return std::real(detail::wright_omega(z, err, res, condn));
      }
  }

/**
 * Return the -1-branch of the Lambert W function of real argument.
 * This is denoted Wm in the DLMF.
 */
template<typename Tp>
  Tp
  lambert_wm(Tp x)
  {
    using Cmplx = std::complex<Tp>;
    const auto s_2pi = emsr::tau_v<Tp>;
    const auto s_i = Cmplx(Tp{0}, Tp{1});
    const auto s_1de = emsr::inv_e_v<Tp>;
    if (x < -s_1de || x > Tp{0})
      return std::numeric_limits<Tp>::quiet_NaN();
    else if (x == -s_1de)
      return Tp{-1};
    else if (x == Tp{0})
      return -std::numeric_limits<Tp>::infinity();
    else
      {
	Cmplx z = std::log(x + s_i * std::numeric_limits<Tp>::epsilon()) - s_i * s_2pi;
	Cmplx err, res, condn;
	return std::real(detail::wright_omega(z, err, res, condn));
      }
  }

} // namespace emsr

#endif // SF_ELEMENTARY_H
