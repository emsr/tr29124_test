// Special functions -*- C++ -*-

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

/** @file emsr/sf_parab_cyl.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// References:
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 9, pp. 355-434, Section 10 pp. 435-478
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 240-245

#ifndef SF_PARAB_CYL_TCC
#define SF_PARAB_CYL_TCC 1

#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{

  /**
   * Return the prefactors used in parabolic cylinder functions.
   */
  template<typename _Tp>
    std::tuple<_Tp, _Tp, _Tp, _Tp>
    parabolic_cylinder_factor(_Tp a)
    {
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto _S_sqrt_pi = _S_pi / _Tp{2};
      auto __2e14p = std::pow(_Tp{2}, 0.25L + 0.5L * a);
      auto __2e34p = std::pow(_Tp{2}, 0.75L + 0.5L * a);
      auto __2e14m = std::pow(_Tp{2}, 0.25L - 0.5L * a);
      auto gamma34p = gamma(0.75L + 0.5L * a);
      auto gamma34m = gamma(0.75L - 0.5L * a);
      auto gamma14p = gamma(0.25L + 0.5L * a);
      auto gamma14m = gamma(0.25L - 0.5L * a);
      auto U0 = _S_sqrt_pi / (__2e14p * gamma34p);
      auto V0 = -_S_sqrt_pi * __2e14m / gamma14p;
      auto Up0 = _S_pi * __2e14p / (gamma34m * gamma34m * gamma14p);
      auto Vp0 = _S_pi * __2e34p / (gamma14m * gamma14m * gamma34p);
      return std::make_tuple(U0, V0, Up0, Vp0);
    }

  /**
   * Return the parabolic cylinder functions by series solution.
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    parabolic_cylinder_series(_Tp a, _Tp z)
    {
      const auto _S_eps = emsr::epsilon(std::real(z));
      constexpr auto _S_max_iter = 1000;
      const auto zz = z * z;
      const auto ezz4 = std::exp(-zz / _Tp{4});
      auto term1 = _Tp{1};
      auto sum1 = term1;
      auto term2 = z;
      auto sum2 = term2;
      for (int k = 1; k < _S_max_iter; ++k)
	{
	  term1 *= (a + _Tp(4 * k - 3) / _Tp{2})
		   * zz / _Tp(2 * k * (2 * k - 1));
	  sum1 += term1;
	  term2 *= (a + _Tp(4 * k - 1) / _Tp{2})
		   * zz / _Tp(2 * k * (2 * k + 1));
	  sum2 += term2;
	  if (std::abs(term1) < _S_eps * std::abs(sum1)
	   && std::abs(term2) < _S_eps * std::abs(sum2))
	    break;
	}
      sum1 *= ezz4;
      sum2 *= ezz4;
      auto fact = parabolic_cylinder_factor(a);
      auto _U = std::get<0>(fact) * sum1 + std::get<2>(fact) * sum2;
      auto _V = std::get<1>(fact) * sum1 + std::get<3>(fact) * sum2;

      return std::make_pair(_U, _V);
    }

  /**
   * Return the parabolic cylinder functions by asymptotic series solution.
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    parabolic_cylinder_asymp(_Tp a, _Tp z)
    {
      const auto _S_eps = emsr::epsilon(std::real(z));
      constexpr auto _S_max_iter = 1000;
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      const auto _S_sqrt_pi = emsr::sqrtpi_v<_Tp>;
      const auto _S_2dsqrt_pi = _Tp{2} / _S_sqrt_pi;
      const auto zz = z * z;
      const auto i2zz = _Tp{1} / (_Tp{2} * z * z);
      const auto ezz4 = std::exp(-zz / _Tp{4});
      const auto pow = std::pow(z, -a - _S_1d2);
      auto term1 = _Tp{1};
      auto sum1 = term1;
      auto term2 = _Tp{1};
      auto sum2 = term2;
      for (auto s = 1; s < _S_max_iter; ++s)
	{
	  term1 *= -(a + _S_1d2 + 2 * (s - 1))
		  * (a + _S_1d2 + 2 * (s - 1) + 1) * i2zz / s;
	  sum1 += term1;
	  term2 *= (a - _S_1d2 + 2 * (s - 1))
		  * (a - _S_1d2 + 2 * (s - 1) + 1) * i2zz / s;
	  sum2 += term2;
	  if (std::abs(term1) < _S_eps * std::abs(sum1)
	   && std::abs(term2) < _S_eps * std::abs(sum2))
	    break;
	}
      auto _U = ezz4 * pow * sum1;
      auto _V = _S_2dsqrt_pi * sum2 / ezz4 / pow / z;

      return std::make_pair(_U, _V);
    }

  /**
   * 
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    parabolic_cylinder(_Tp a, _Tp z)
    {
      constexpr auto _S_magic_switch = _Tp{10};
      if (std::abs(z) < _S_magic_switch)
	return parabolic_cylinder_series(a, z);
      else
	return parabolic_cylinder_asymp(a, z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    parabolic_cyl_u(_Tp a, _Tp z)
    {
      if (std::isnan(a) || std::isnan(z))
	return emsr::quiet_NaN<_Tp>();
      else
        return parabolic_cylinder(a, z).first;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    parabolic_cyl_v(_Tp a, _Tp z)
    {
      if (std::isnan(a) || std::isnan(z))
	return emsr::quiet_NaN<_Tp>();
      else
        return parabolic_cylinder(a, z).second;
    }

  /**
   *  
   */
  template<typename _Tp>
    _Tp
    parabolic_cyl_w(_Tp a, _Tp z)
    {
      return _Tp{0};
    }

} // namespace detail
} // namespace emsr

#endif // SF_PARAB_CYL_TCC
