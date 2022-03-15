
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

/** @file bits/sf_trig.tcc
 * This is an internal header file, included by other library headers.
 * You should not attempt to use it directly.
 */

#ifndef SF_TRIG_TCC
#define SF_TRIG_TCC 1

#include <complex>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{
  /**
   * A type describing a cosine and a sine value.
   * A return for sincos-type functions.
   */
  template<typename Tp>
    struct sincos_t
    {
      Tp sin_v;
      Tp cos_v;
    };

namespace detail
{

  /**
   * Return the reperiodized sine of argument x:
   * @f[
   *   \mathrm{sin_\pi}(x) = \sin(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    sin_pi(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Val>;
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return -sin_pi(-x);
      else if (x < Tp{0.5L})
	return std::sin(x * s_pi);
      else if (x < Tp{1})
	return std::sin((Tp{1} - x) * s_pi);
      else
	{
	  auto nu = std::floor(x);
	  auto arg = x - nu;
	  auto sign = (int(nu) & 1) == 1 ? -1 : +1;
	  auto sinval = (arg < Tp{0.5L})
			? sin_pi(arg)
			: sin_pi(Tp{1} - arg);
	  return sign * sinval;
	}
    }

  /**
   * Return the reperiodized hyperbolic sine of argument x:
   * @f[
   *   \mathrm{sinh_\pi}(x) = \sinh(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    sinh_pi(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Val>;
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return -sinh_pi(-x);
      else
	return std::sinh(s_pi * x);
    }

  /**
   * Return the reperiodized cosine of argument x:
   * @f[
   *   \mathrm{cos_\pi}(x) = \cos(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    cos_pi(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Val>;
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return cos_pi(-x);
      else if (x < Tp{0.5L})
	return std::cos(x * s_pi);
      else if (x < Tp{1})
	return -std::cos((Tp{1} - x) * s_pi);
      else
	{
	  auto nu = std::floor(x);
	  auto arg = x - nu;
	  auto sign = (int(nu) & 1) == 1 ? -1 : +1;
	  return sign * cos_pi(arg);
	}
    }

  /**
   * Return the reperiodized hyperbolic cosine of argument x:
   * @f[
   *   \mathrm{cosh_\pi}(x) = \cosh(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    cosh_pi(Tp x)
    {
      using Val = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Val>;
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return cosh_pi(-x);
      else
	return std::cosh(s_pi * x);
    }

  /**
   * Return the reperiodized tangent of argument x:
   * @f[
   *   \tan_pi(x) = \tan(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    tan_pi(Tp x)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      return std::tan(s_pi * (x - std::floor(x)));
    }

  /**
   * Return the reperiodized hyperbolic tangent of argument x:
   * @f[
   *   \mathrm{tanh_\pi}(x) = \tanh(\pi x)
   * @f]
   */
  template<typename Tp>
    Tp
    tanh_pi(Tp x)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      return std::tanh(s_pi * x);
    }

  /**
   * Return the reperiodized sine of complex argument z:
   * @f[
   *   \mathrm{sin_\pi}(z) = \sin(\pi z)
   *     = \mathrm{sin_\pi}(x) \mathrm{cosh_\pi}(y)
   *   + i \mathrm{cos_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    sin_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      return sin_pi(x) * std::cosh(s_pi * y)
	   + s_i * cos_pi(x) * std::sinh(s_pi * y);
    }

  /**
   * Return the reperiodized hyperbolic sine of complex argument z:
   * @f[
   *   \mathrm{sinh_\pi}(z) = \sinh(\pi z)
   *     = \mathrm{\sinh_\pi}(x) \mathrm{cos_\pi}(y)
   *   + i \mathrm{\cosh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    sinh_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      return std::sinh(s_pi * x) * cos_pi(y)
	+ s_i * std::cosh(s_pi * x) * sin_pi(y);
    }

  /**
   * Return the reperiodized cosine of complex argument z:
   * @f[
   *    \mathrm{cos_\pi}(z) = \cos(\pi z)
   *       = \mathrm{cos_\pi}(x) \mathrm{cosh_\pi}(y)
   *     - i \mathrm{sin_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    cos_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      return cos_pi(x) * std::cosh(s_pi * y)
	   - s_i * sin_pi(x) * std::sinh(s_pi * y);
    }

  /**
   * Return the reperiodized hyperbolic cosine of complex argument z:
   * @f[
   *    \mathrm{cosh_\pi}(z) = \mathrm{cosh_\pi}(z)
   *       = \mathrm{cosh_\pi}(x) \mathrm{cos_\pi}(y)
   *     + i \mathrm{sinh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    cosh_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      return std::cosh(s_pi * x) * cos_pi(y)
	   + s_i * std::sinh(s_pi * x) * sin_pi(y);
    }

  /**
   * Return the reperiodized tangent of complex argument z:
   * @f[
   *   \mathrm{tan_\pi}(z) = \tan(\pi z)
   *     = \frac{\mathrm{tan_\pi}(x) + i \mathrm{tanh_\pi}(y)}
   *            {1 - i \mathrm{tan_\pi}(x) \mathrm{tanh_\pi}(y)}
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    tan_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      auto tan = tan_pi(x);
      auto tanh = std::tanh(s_pi * y);
      return (tan + s_i * tanh) / (Val{1} - s_i * tan * tanh);
    }

  /**
   * Return the reperiodized hyperbolic tangent of complex argument z:
   * @f[
   *   \mathrm{tanh_\pi}(z) = \tanh(\pi z)
   *     = \frac{\mathrm{tanh_\pi}(x) + i \mathrm{tan_\pi}(y)}
   *            {1 + i \mathrm{tanh_\pi}(x) \mathrm{tan_\pi}(y)}
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    tanh_pi(std::complex<Tp> z)
    {
      using Val = Tp;
      using Real = emsr::num_traits_t<Val>;
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_i = std::complex<Tp>{0, 1};
      auto x = std::real(z);
      auto y = std::imag(z);
      auto tanh = std::tanh(s_pi * x);
      auto tan = tan_pi(y);
      return (tanh + s_i * tan) / (Val{1} + s_i * tanh * tan);
    }

  /**
   * 
   */
  template<typename Tp>
    inline emsr::sincos_t<Tp>
    sincos(Tp x)
    { return emsr::sincos_t<Tp>{std::sin(x), std::cos(x)}; }

  /**
   * 
   */
  template<>
    inline emsr::sincos_t<float>
    sincos(float x)
    {
      float sinx, cosx;
      __builtin_sincosf(x, &sinx, &cosx);
      return emsr::sincos_t<float>{sinx, cosx};
    }

  /**
   * 
   */
  template<>
    inline emsr::sincos_t<double>
    sincos(double x)
    {
      double sinx, cosx;
      __builtin_sincos(x, &sinx, &cosx);
      return emsr::sincos_t<double>{sinx, cosx};
    }

  /**
   * 
   */
  template<>
    inline emsr::sincos_t<long double>
    sincos(long double x)
    {
      long double sinx, cosx;
      __builtin_sincosl(x, &sinx, &cosx);
      return emsr::sincos_t<long double>{sinx, cosx};
    }

  /**
   * Reperiodized sincos.
   */
  template<typename Tp>
    emsr::sincos_t<Tp>
    sincos_pi(Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_NaN = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return emsr::sincos_t<Tp>{s_NaN, s_NaN};
      else if (x < Tp{0})
	{
	  emsr::sincos_t<Tp> tempsc = sincos_pi(-x);
	  return emsr::sincos_t<Tp>{-tempsc.sin_v, tempsc.cos_v};
	}
      else if (x < Tp{0.5L})
	return sincos(s_pi * x);
      else if (x < Tp{1})
	{
	  emsr::sincos_t<Tp> tempsc = sincos(s_pi * (Tp{1} - x));
	  return emsr::sincos_t<Tp>{tempsc.sin_v, -tempsc.cos_v};
	}
      else
	{
	  auto nu = std::floor(x);
	  auto arg = x - nu;
	  auto sign = (int(nu) & 1) == 1 ? Tp{-1} : Tp{+1};

	  auto sinval = (arg < Tp{0.5L})
			? std::sin(s_pi * arg)
			: std::sin(s_pi * (Tp{1} - arg));
	  auto cosval = std::cos(s_pi * arg);
	  return emsr::sincos_t<Tp>{sign * sinval,
					    sign * cosval};
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_TRIG_TCC
