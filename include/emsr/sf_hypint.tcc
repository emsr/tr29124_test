
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

/** @file bits/sf_hypint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SF_HYPINT_TCC
#define SF_HYPINT_TCC 1

#include <stdexcept>
#include <complex>

namespace emsr
{
namespace detail
{

  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by continued fraction for positive argument.
   */
  ////FIXME!!!!
  template<typename _Tp>
    void
    chshint_cont_frac(_Tp t, _Tp& _Chi, _Tp& _Shi)
    {
      const unsigned int s_max_iter = 100;
      const auto s_eps = _Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};

      // Evaluate Chi and Shi by Lentz's modified method of continued fracions.
      std::complex<_Tp> b(_Tp{1}, t);
      std::complex<_Tp> c(_Tp{1} / s_fp_min);
      std::complex<_Tp> d(_Tp{1} / b);
      std::complex<_Tp> h(d);
      unsigned int i = 2;
      while (true)
	{
	  auto a = -_Tp(i - 1) * _Tp(i - 1);
	  b += _Tp{2};
	  d = _Tp{1} / (a * d + b);
	  c = b + a / c;
	  auto del = c * d;
	  h *= del;
	  if (std::abs(del.real() - _Tp{1}) + std::abs(del.imag()) < s_eps)
	    break;
	  if (i > s_max_iter)
	    throw std::runtime_error("chshint_cont_frac: Continued fraction evaluation failed");
	  ++i;
	}
      h *= std::polar(_Tp{1}, -t);
      _Chi = -h.real();
      _Shi = s_pi_2 + h.imag();

      return;
    }


  /**
   *  @brief This function computes the hyperbolic cosine @f$ Chi(x) @f$
   *    and hyperbolic sine @f$ Shi(x) @f$ integrals
   *    by series summation for positive argument.
   */
  template<typename _Tp>
    void
    chshint_series(_Tp t, _Tp& _Chi, _Tp& _Shi)
    {
      const auto s_max_iter = 100;
      const auto s_eps = _Tp{5} * emsr::epsilon(t);
      const auto s_fp_min = emsr::lim_min(t);
      const auto s_gamma_e = emsr::egamma_v<_Tp>;

      // Evaluate Chi and Shi by series simultaneously.
      _Tp _Csum(0), _Ssum(0);
      if (t * t < s_fp_min)
	{
	  // Avoid underflow.
	  _Csum = _Tp{0};
	  _Ssum = t;
	}
      else
	{
	  // Evaluate Shi and Chi by series expansion.
	  _Tp sum(0);
	  _Tp fact(1);
	  auto odd = true;
	  auto k = 1;
	  while (true)
	    {
	      fact *= t / k;
	      _Tp term = fact / k;
	      sum += term;
	      _Tp err = term / std::abs(sum);
	      if (odd)
		{
		  _Ssum = sum;
		  sum = _Csum;
		}
	      else
		{
		  _Csum = sum;
		  sum = _Ssum;
		}
	      if (err < s_eps)
		break;
	      odd = !odd;
	      ++k;
	      if (k > s_max_iter)
		throw std::runtime_error("chshint_series: Series evaluation failed");
	    }
	}
      _Chi = s_gamma_e + std::log(t) + _Csum;
      _Shi = _Ssum;

      return;
    }


  /**
   *  @brief This function returns the hyperbolic cosine @f$ Ci(x) @f$
   *    and hyperbolic sine @f$ Si(x) @f$ integrals as a pair.
   *
   *  The hyperbolic cosine integral is defined by:
   *  @f[
   *      Chi(x) = \gamma_E + \log(x) + \int_0^x dt \frac{\cosh(t) - 1}{t}
   *  @f]
   *
   *  The hyperbolic sine integral is defined by:
   *  @f[
   *      Shi(x) = \int_0^x dt \frac{\sinh(t)}{t}
   *  @f]
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    chshint(_Tp x, _Tp& _Chi, _Tp& _Shi)
    {
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return std::make_pair(s_NaN, s_NaN);

      auto t = std::abs(x);
      if (t == _Tp{0})
	{
	  _Chi = -emsr::infinity(x);
	  _Shi = _Tp{0};
	}
      else if (t > _Tp{2})
	chshint_cont_frac(t, _Chi, _Shi);
      else
	chshint_series(t, _Chi, _Shi);

      if (x < _Tp{0})
	_Shi = -_Shi;

      return std::make_pair(_Chi, _Shi);
    }

} // namespace detail
} // namespace emsr

#endif // SF_HYPINT_TCC
