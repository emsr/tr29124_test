
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_trigint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_TRIGINT_TCC
#define SF_TRIGINT_TCC 1

#include <stdexcept>
#include <complex>

namespace emsr
{
namespace detail
{

  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by continued fraction
   *         for positive argument.
   */
  template<typename _Tp>
    void
    sincosint_cont_frac(_Tp t, _Tp& _Si, _Tp& _Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = _Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};

      // Evaluate Ci and Si by Lentz's modified method of continued fractions.
      std::complex<_Tp> b(_Tp{1}, t);
      std::complex<_Tp> c(_Tp{1} / s_fp_min);
      std::complex<_Tp> d(_Tp{1} / b);
      std::complex<_Tp> h(d);
      int i = 2;
      while (true)
	{
	  auto a = -_Tp(i - 1) * _Tp(i - 1);
	  b += _Tp{2};
	  d = _Tp{1} / (a * d + b);
	  c = b + a / c;
	  std::complex<_Tp> del = c * d;
	  h *= del;
	  if (std::abs(del - _Tp{1}) < s_eps)
	    break;
	  if (i > s_max_iter)
	    throw std::runtime_error("sincosint_cont_frac: continued fraction evaluation failed");
	  ++i;
	}
      h *= std::polar(_Tp{1}, -t);
      _Ci = -h.real();
      _Si = s_pi_2 + h.imag();

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by series summation
   *         for positive argument.
   */
  template<typename _Tp>
    void
    sincosint_series(_Tp t, _Tp& _Si, _Tp& _Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = _Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_gamma_e = emsr::egamma_v<_Tp>;

      // Evaluate Ci and Si by series simultaneously.
      _Tp sumc(0), sums(0);
      if (t * t < s_fp_min)
	{
	  // Avoid underflow.
	  sumc = _Tp{0};
	  sums = t;
	}
      else
	{
	  // Evaluate Si and Ci by series expansion.
	  _Tp sum(0);
	  _Tp sign(1), fact(1);
	  bool odd = true;
	  unsigned int k = 1;
	  while (true)
	    {
	      fact *= t / k;
	      _Tp term = fact / k;
	      sum += sign * term;
	      _Tp err = term / std::abs(sum);
	      if (odd)
		{
		  sign = -sign;
		  sums = sum;
		  sum = sumc;
		}
	      else
		{
		  sumc = sum;
		  sum = sums;
		}
	      if (err < s_eps)
		break;
	      odd = !odd;
	      ++k;
	      if (k > s_max_iter)
		throw std::runtime_error("sincosint_series: series evaluation failed");
	    }
	}
      _Si = sums;
      _Ci = s_gamma_e + std::log(t) + sumc;

      return;
    }


  /**
   *  @brief This function computes the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals by asymptotic series summation
   *         for positive argument.
   *
   *   The asymptotic series is very good for x > 50.
   */
  template<typename _Tp>
    void
    sincosint_asymp(_Tp t, _Tp& _Si, _Tp& _Ci)
    {
      const auto s_max_iter = 100;
      const auto s_eps = _Tp{5} * emsr::epsilon(t);
      const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};

      auto invt = _Tp{1} / t;
      auto term = _Tp{1}; // 0!
      auto sume = _Tp{term};
      term *= invt; // 1! / t
      auto sumo = _Tp{term};
      auto sign = _Tp{1};
      auto even = true;
      auto k = 2;
      while (true)
	{
	  term *= k * invt;

	  if (even)
	    {
	      sign = -sign;
	      sume += sign * term;

	    }
	  else
	    {
	      sumo += sign * term;
	      if (term / std::abs(sumo) < s_eps)
		break;
	    }

	  even = !even;

	  if (k > s_max_iter)
	    throw std::runtime_error("sincosint_asymp: Series evaluation failed");
	  ++k;
	}

      auto sint = std::sin(t);
      auto cost = std::cos(t);
      _Si = s_pi_2
	  - cost * invt * sume
	  - sint * invt * sumo;
      _Ci = sint * invt * sume
	  - cost * invt * sumo;

      return;
    }


  /**
   *  @brief This function returns the sine @f$ Si(x) @f$
   *         and cosine @f$ Ci(x) @f$ integrals as a @c pair.
   *
   *  The sine integral is defined by:
   *  @f[
   *      Si(x) = \int_0^x dt \frac{\sin(t)}{t}
   *  @f]
   *
   *  The cosine integral is defined by:
   *  @f[
   *      Ci(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cos(t) - 1}{t}
   *  @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    sincosint(_Tp x)
    {
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return std::make_pair(s_NaN, s_NaN);

      auto t = std::abs(x);
      _Tp _Ci, _Si;
      if (t == _Tp{0})
	{
	  _Si = _Tp{0};
	  _Ci = -emsr::infinity(x);
	}
      else if (t > _Tp{1000}) // Check this!
	sincosint_asymp(t, _Si, _Ci);
      else if (t > _Tp{2})
	sincosint_cont_frac(t, _Si, _Ci);
      else
	sincosint_series(t, _Si, _Ci);

      if (x < _Tp{0})
	_Si = -_Si;

      return std::make_pair(_Si, _Ci);
    }

} // namespace detail
} // namespace emsr

#endif // SF_TRIGINT_TCC
