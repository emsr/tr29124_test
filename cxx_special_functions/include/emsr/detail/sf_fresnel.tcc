
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

/** @file bits/sf_fresnel.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_FRESNEL_TCC
#define SF_FRESNEL_TCC 1

#include <stdexcept>
#include <complex>

namespace emsr
{
namespace detail
{

  /**
   *  @brief This function returns the Fresnel cosine and sine integrals
   *    as a pair by series expansion for positive argument.
   */
  template <typename Tp>
    void
    fresnel_series(const Tp ax, Tp & _Cf, Tp & _Sf)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(ax);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = s_pi / Tp{2};

      // Evaluate S and C by series expansion.
      auto sum = Tp{0};
      auto _Ssum = Tp{0};
      auto _Csum = ax;
      auto sign = Tp{1};
      auto fact = s_pi_2 * ax * ax;
      auto odd = true;
      auto term = ax;
      auto n = 3;
      auto k = 0;
      for (k = 1; k <= s_max_iter; ++k)
	{
	  term *= fact / k;
	  sum += sign * term / n;
	  Tp test = std::abs(sum) * s_eps;
	  if (odd)
	    {
	      sign = -sign;
	      _Ssum = sum;
	      sum = _Csum;
	    }
	  else
	    {
	      _Csum = sum;
	      sum = _Ssum;
	    }

	  if (term < test)
	    break;

	  odd = ! odd;

	  n += 2;
	}
      if (k > s_max_iter)
	throw std::runtime_error("fresnel_series: series evaluation failed");

      _Cf = _Csum;
      _Sf = _Ssum;

      return;
    }

  /**
   *  @brief This function computes the Fresnel cosine and sine integrals
   *    by continued fractions for positive argument.
   */
  template <typename Tp>
    void
    fresnel_cont_frac(const Tp ax, Tp & _Cf, Tp & _Sf)
    {
      const auto s_max_iter = 100;
      const auto s_eps = Tp{5} * emsr::epsilon(ax);
      const auto s_fp_min = emsr::lim_min(ax);
      const auto s_pi = emsr::pi_v<Tp>;

      // Evaluate S and C by Lentz's complex continued fraction method.
      const auto pix2 = s_pi * ax * ax;
      std::complex<Tp> b(Tp{1}, -pix2);
      std::complex<Tp> cc(Tp{1} / s_fp_min, Tp{0});
      auto h = Tp{1} / b;
      auto d = h;
      auto n = -1;
      auto k = 0;
      for (k = 2; k <= s_max_iter; ++k)
	{
	  n += 2;
	  const auto a = -Tp(n * (n + 1));
	  b += Tp{4};
	  d = Tp{1} / (a * d + b);
	  cc = b + a / cc;
	  const auto del = cc * d;
	  h *= del;
	  if (std::abs(del.real() - Tp{1})
	    + std::abs(del.imag()) < s_eps)
	    break;
	}
      if (k > s_max_iter)
	throw std::runtime_error("fresnel_cont_frac: continued fraction evaluation failed");

      h *= std::complex<Tp>(ax, -ax);
      auto phase = std::polar(Tp{1}, pix2/Tp{2});
      auto cs = std::complex<Tp>(Tp{0.5L}, Tp{0.5L})
		* (Tp{1} - phase * h);
      _Cf = cs.real();
      _Sf = cs.imag();

      return;
    }


  /**
   * @brief Return the Fresnel cosine and sine integrals
   * as a complex number $f[ C(x) + iS(x) $f].
   *
   * The Fresnel cosine integral is defined by:
   * @f[
   * 	 C(x) = \int_0^x \cos(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * The Fresnel sine integral is defined by:
   * @f[
   * 	 S(x) = \int_0^x \sin(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param x The argument
   */
  template <typename Tp>
    std::complex<Tp>
    fresnel(const Tp x)
    {
      const auto s_fp_min = emsr::lim_min(x);
      const auto s_x_min = Tp{1.5L};
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return std::complex<Tp>{s_NaN, s_NaN};

      auto _Cf = Tp{0};
      auto _Sf = Tp{0};

      const Tp ax = std::abs(x);
      if (ax < std::sqrt(s_fp_min))
	{
	  _Cf = ax;
	  _Sf = Tp{0};
	}
      else if (ax < s_x_min)
	fresnel_series(ax, _Cf, _Sf);
      else
	fresnel_cont_frac(ax, _Cf, _Sf);

      if (x < Tp{0})
	{
	  _Cf = -_Cf;
	  _Sf = -_Sf;
	}

      return std::complex<Tp>(_Cf, _Sf);
    }

} // namespace detail
} // namespace emsr

#endif // SF_FRESNEL_TCC
