
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

/** @file bits/sf_elementary.tcc
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef SF_ELEMENTARY_TCC
#define SF_ELEMENTARY_TCC 1

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

template<typename Tp>
  std::complex<Tp>
  wright_omega(const std::complex<Tp>& z,
		 std::complex<Tp>& err, std::complex<Tp>& res,
		 std::complex<Tp>& condn)
  {
    using Cmplx = std::complex<Tp>;
    const auto s_NaN = emsr::quiet_NaN(z.real());
    const auto s_eps = emsr::epsilon(z.real());
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_i = Cmplx(Tp{0}, Tp{1});
    const auto s_0 = Cmplx{};
    auto [x, y] = reinterpret_cast<const Tp(&)[2]>(z);
    auto ympi = y - s_pi;
    auto yppi = y + s_pi;
    const auto s_near = Tp{0.1};
    err = s_0;
    res = s_0;

    if (std::isnan(x) || std::isnan(y))
      return Cmplx(s_NaN, s_NaN);
    else if (std::isinf(x) && (x < Tp{0})
	  && (-s_pi < y) && (y <= s_pi))
      { // Signed zeros between branches.
	Cmplx w;
	if (std::abs(y) <= s_pi / Tp{2})
	  w = +Tp{0};
	else
	  w = -Tp{0};

	if (y < Tp{0})
	  w += -s_i * Tp{0};

	return w;
      }
    else if (std::isinf(x) || std::isinf(y))
      return Cmplx(x, y);
    else if (x == Tp{-1} && std::abs(y) == s_pi)
      return Cmplx(Tp{-1}, Tp{0});
    else
      {
	//  Choose approximation based on region.

	Cmplx w;
	if ((Tp{-2} < x && x <= Tp{1}
	  && Tp{1} < y && y < Tp{2} * s_pi))
	  {
	    // Region 1: upper branch point.
	    // Series about z = -1 + i*pi.
	    const auto dz = z + Tp{1} - s_i * s_pi;
	    const auto pz = std::conj(std::sqrt(std::conj(Tp{2} * dz)));

	    w = Tp{-1}
		+ (s_i
		+ (Tp{1} / Tp{3}
		+ (Tp{-1} / Tp{36} * s_i
		+ (Tp{1} / Tp{270} + Tp{1} / Tp{4320} * s_i * pz)
		* pz) * pz) * pz) * pz;
	  }
	else if ((Tp{-2} < x && x <= Tp{1}
	       && Tp{-2} * s_pi < y && y < Tp{-1}))
	  {
	    // Region 2: lower branch point.
	    // Series about z = -1 - i*pi.
	    const auto dz = z + Tp{1} + s_i * s_pi;
	    const auto pz = std::conj(std::sqrt(std::conj(Tp{2} * dz)));

	    w = Tp{-1}
		+ (-s_i
		+ (Tp{1} / Tp{3}
		+ (Tp{1} / Tp{36} * s_i
		+ (Tp{1} / Tp{270} - Tp{1} / Tp{4320} * s_i * pz)
		* pz) * pz) * pz) * pz;
	  }
	else if (x <= Tp{-2} && -s_pi < y && y <= s_pi)
	  {
	    // Region 3: between branch cuts.
	    // Series: About -infinity.
	    const auto pz = std::exp(z);
	    w = (Tp{1}
		+ (Tp{-1}
		+ (Tp{3} / Tp{2}
		+ (Tp{-8} / Tp{3}
		+ Tp{125} / Tp{24} * pz) * pz) * pz) * pz) * pz;
	  }
	else if (((Tp{-2} < x) && (x <= Tp{1})
	       && (Tp{-1} <= y) && (y <= Tp{1}))
		|| ((Tp{-2} < x)
		 && (x - Tp{1}) * (x - Tp{1}) + y * y
			 <= s_pi * s_pi))
	  {
	    // Region 4: Mushroom.
	    // Series about z = 1.
	    const auto pz = z - Tp{1};
	    w = Tp{1} / Tp{2} + Tp{1} / Tp{2} * z
		+ (Tp{1} / Tp{16}
		+ (Tp{-1} / Tp{192}
		+ (Tp{-1} / Tp{3072}
		+ Tp{13} / Tp{61440} * pz) * pz) * pz) * pz * pz;
	  }
	else if (x <= -Tp{3} / Tp{2}
		 && s_pi < y
		 && y - s_pi <= Tp{-3} / Tp{4} * (x + Tp{1}))
	  {
	    // Region 5: Top wing.
	    // Negative log series.
	    const auto t = z - s_i * s_pi;
	    const auto pz = std::log(-t);
	    w = ((Tp{1}
		  + (-Tp{3} / Tp{2}
		  + Tp{1} / Tp{3} * pz) * pz) * pz
		 + ((Tp{-1}
		  + Tp{1} / Tp{2} * pz) * pz
		  + (pz + (-pz + t) * t) * t) * t)
		/ (t * t * t);
	  }
	else if (x <= -Tp{3} / Tp{2}
		 && Tp{3} / Tp{4} * (x + Tp{1}) < y + s_pi
				    && y + s_pi <= Tp{0})
	  {
	    // Region 6: Bottom wing.
	    // Negative log series.
	    const auto t = z + s_i * s_pi;
	    const auto pz = std::log(-t);
	    w = ((Tp{1}
		 + (Tp{-3} / Tp{2}
		 + Tp{1} / Tp{3} * pz) * pz) * pz
		+ ((Tp{-1}
		 + Tp{1} / Tp{2} * pz) * pz
		 + (pz + (-pz + t) * t) * t) * t)
		/ (t * t * t);
	  }
	else
	  {
	    // Region 7: Everywhere else.
	    // Series solution about infinity.
	    const auto pz = std::log(z);
	    w = ((Tp{1}
		 + (Tp{-3} / Tp{2}
		 + Tp{1} / Tp{3} * pz) * pz) * pz
		+ ((Tp{-1}
		 + Tp{1} / Tp{2} * pz) * pz
		 + (pz + (-pz + z) * z) * z) * z)
		/ (z * z * z);
	  }

	auto sgn = Tp{0};
	auto zr = z;
	if (x <= Tp{-1} + s_near
	    && (std::abs(ympi) <= s_near || std::abs(yppi) <= s_near))
	  {
	    sgn = Tp{-1};
	    // Regularize if near branch cuts.
	    if (std::abs(ympi) <= s_near)
	      {
		// Recompute ympi with directed rounding.
		fesetround(FE_UPWARD);

		ympi = y - s_pi;

		if (ympi <= Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    ympi = y - s_pi;
		  }

		zr = Cmplx(x, ympi);

		// Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	    else
	      {
		// Recompute yppi with directed rounding.
		fesetround(FE_UPWARD);

		yppi = y + s_pi;

		if (yppi <= Tp{0})
		  {
        	    fesetround(FE_DOWNWARD);
        	    yppi = y + s_pi;
		  }

		zr = Cmplx(x, yppi);

		//  Return rounding to default.
		fesetround(FE_TONEAREST);
	      }
	  }
	else
	  sgn = Tp{+1};

	w *= sgn;

	while (true)
	  {
	    const auto res = zr - sgn * w - std::log(w);
	    const auto wp1 = sgn * w + Tp{1};
	    const auto yy = Tp{2} * wp1 * (wp1 + Tp{2} / Tp{3} * res);
	    const auto err = res / wp1 * (yy - res)
				/ (yy - Tp{2} * res);
	    w *= Tp{1} + err;
	    const auto res4 = std::pow(std::abs(res), Tp{4});
	    const auto wpol = Tp{-1} + w * (Tp{-8} + Tp{2} * w);
	    const auto test = std::abs(res4 * wpol);
	    if (test < s_eps * Tp{72} * std::pow(std::abs(wp1), Tp{6}))
	      break;
	  }

	// Undo regularization.
	w *= sgn;

	// Provide condition number estimate.
	condn = zr / (Tp{1} + w);

	return w;
      }
  }

} // namespace detail
} // namespace emsr

#endif // SF_ELEMENTARY_TCC
