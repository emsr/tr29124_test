
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

/** @file bits/sf_bessel.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
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

#ifndef SF_BESSEL_TCC
#define SF_BESSEL_TCC 1

#include <stdexcept>
#include <complex>
#include <utility> // For exchange.

#include <emsr/fp_type_util.h>
#include <emsr/sf_trig.h> // sin_pi, cos_pi, polar_pi
#include <emsr/specfun_state.h>

namespace emsr
{
namespace detail
{

  /**
   * @brief This routine returns the cylindrical Bessel functions
   * 	    of order @f$ \nu @f$: @f$ J_{\nu}(z) @f$ or @f$ I_{\nu}(z) @f$
   * 	    by series expansion.
   *
   * The modified cylindrical Bessel function is:
   * @f[
   *  Z_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{\sigma^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   * where @f$ \sigma = +1 @f$ or @f$ -1 @f$ for
   * @f$ Z = I @f$ or @f$ J @f$ respectively.
   *
   * See Abramowitz & Stegun, 9.1.10
   * 	 Abramowitz & Stegun, 9.6.7
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Equation 9.1.10 p. 360 and Equation 9.6.10 p. 375
   *
   * @param  nu  The order of the Bessel function.
   * @param  x   The argument of the Bessel function.
   * @param  sgn  The sign of the alternate terms
   * 		    -1 for the Bessel function of the first kind.
   * 		    +1 for the modified Bessel function of the first kind.
   * @param  max_iter  The maximum number of iterations for sum.
   * @return  The output Bessel function.
   */
  template<typename _Tnu, typename _Tp>
    constexpr _Tp
    cyl_bessel_ij_series(_Tnu nu, _Tp x, int sgn,
			   unsigned int max_iter)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, _Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon<_Real>();
      if (std::abs(x) < s_eps)
	{
	  if (nu == _Tnu{0})
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	{
	  const auto x2 = x / _Real{2};

	  const auto xx4 = _Val(sgn) * x2 * x2;
	  auto _Jn = _Val{1};
	  auto term = _Val{1};
	  for (unsigned int i = 1; i < max_iter; ++i)
	    {
	      term *= xx4 / (_Val(i) * (_Val(nu) + _Val(i)));
	      _Jn += term;
	      if (std::abs(term / _Jn) < s_eps)
		break;
	    }

	  auto fact = _Val(nu) * std::log(x2);
	  fact -= log_gamma(_Real{1} + nu);
	  fact = std::exp(fact);
	  return fact * _Jn;
	}
    }

  /**
   * A type for Bessel asymptotic sums.
   */
  template<typename _Tnu, typename _Tp>
    struct cyl_bessel_asymp_sums_t
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, _Tp>;
      _Val _Psum;
      _Val _Qsum;
      _Val _Rsum;
      _Val _Ssum;
    };

  /**
   * @brief This routine computes the asymptotic cylindrical Bessel
   * 	    and Neumann functions of order nu: @f$ J_{\nu}(z) @f$,
   * 	    @f$ N_{\nu}(z) @f$.  Use this for @f$ z >> \nu^2 + 1 @f$.
   *
   * @f[
   *   J_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  - \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * and
   * @f[
   *   N_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  + \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * where @f$ \omega = z - \nu\pi/2 - \pi/4 @f$ and
   * @f[
   *   a_{k}(\nu) = \frac{(4\nu^2 - 1^2)(4\nu^2 - 3^2)...(4\nu^2 - (2k-1)^2)}
   *                     {8^k k!}
   * @f]
   * The derivatives are also computed:
   * @f[
   *   J'_{\nu}(z) = -\left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k}(\nu)}{z^{2k}}
   *  + \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * and
   * @f[
   *   N'_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k}(\nu)}{z^{2k}}
   *  - \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{b_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * where
   * @f[
   *   b_{k}(\nu) = \frac{(4\nu^2 - 1^2)(4\nu^2 - 3^2)...(4\nu^2 - (2k-3)^2)}
   *                     {8^k k!}
   * @f]
   *
   * These sums work everywhere but on the negative real axis:
   * @f$ |ph(z)| < \pi - \delta @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *  (2) Digital Library of Mathematical Fundtions
   *      https://dlmf.nist.gov/10.17
   *
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename _Tp>
    constexpr cyl_bessel_asymp_sums_t<_Tnu, _Tp>
    cyl_bessel_asymp_sums(_Tnu nu, _Tp x, int sgn)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, _Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using bess_t = cyl_bessel_asymp_sums_t<_Tnu, _Tp>;
      const auto s_eps = emsr::epsilon<_Real>();
      const auto __2nu = _Real{2} * nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto r8x = _Tp{1} / (_Real{8} * x);
      const auto nu_min = std::real(nu / _Real{2});
      const auto nu_max = std::abs(_Real{100} * (nu + _Tnu{1}));
      auto k = 0;
      auto bk_xk = _Val{1};
      auto _Rsum = bk_xk;
      auto ak_xk = _Val{1};
      auto _Psum = ak_xk;
      auto convP = false;
      ++k;
      auto __2km1 = 1;
      bk_xk *= (__4nu2 + 3) * r8x;
      auto _Ssum = bk_xk;
      ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * r8x;
      auto _Qsum = ak_xk;
      auto convQ = false;
      auto ak_xk_prev = std::abs(ak_xk);
      do
	{
	  ++k;
	  auto rk8x = r8x / _Real(k);
	  __2km1 += 2;
	  bk_xk = sgn * (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk * rk8x;
	  _Rsum += bk_xk;
	  ak_xk *= sgn * (__2nu - __2km1) * (__2nu + __2km1) * rk8x;
	  if (k > nu_min && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Psum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convP = std::abs(ak_xk) < s_eps * std::abs(_Psum);

	  ++k;
	  rk8x = r8x / _Real(k);
	  __2km1 += 2;
	  bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk * rk8x;
	  _Ssum += bk_xk;
	  ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * rk8x;
	  if (k > nu_min && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Qsum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convQ = std::abs(ak_xk) < s_eps * std::abs(_Qsum);

	  if (convP && convQ)
	    break;
	}
      while (k < nu_max);

      return bess_t{_Psum, _Qsum, _Rsum, _Ssum};
    }

  /**
   *
   */
  template<typename _Tnu, typename _Tp>
    constexpr emsr::cyl_bessel_t<_Tnu, _Tp, _Tp>
    cyl_bessel_jn_asymp(_Tnu nu, _Tp x)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, _Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using bess_t = emsr::cyl_bessel_t<_Tnu, _Tp, _Tp>;
      const auto s_pi = emsr::pi_v<_Real>;
      const auto s_pi_2 = emsr::pi_v<_Real> / _Real{2};

      const auto sums = cyl_bessel_asymp_sums(nu, x, -1);

      const auto omega = x - (nu + _Real{0.5L}) * s_pi_2;
      const auto c = std::cos(omega);
      const auto s = std::sin(omega);

      const auto coef = std::sqrt(_Real{2} / (s_pi * x));
      return bess_t{nu, x,
		 coef * (c * sums._Psum - s * sums._Qsum),
		-coef * (s * sums._Rsum + c * sums._Ssum),
		 coef * (s * sums._Psum + c * sums._Qsum),
		 coef * (c * sums._Rsum - s * sums._Ssum)};
    }

  /**
   * @brief Compute the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 	 \left[\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   *     \left[\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are exquisite.
   *
   * @param mu The input parameter of the gamma functions.
   * @return  An output structure containing four gamma functions.
   */
  template<typename _Tp>
    emsr::gamma_temme_t<_Tp>
    gamma_temme(_Tp mu)
    {
      using gammat_t = emsr::gamma_temme_t<_Tp>;
      const auto s_eps = emsr::epsilon(mu);
      const auto s_gamma_E = emsr::egamma_v<_Tp>;

      if (std::abs(mu) < s_eps)
	return gammat_t{mu, _Tp{1}, _Tp{1}, -s_gamma_E, _Tp{1}};
      else
	{
	  _Tp gamp, gamm;
	  if (std::real(mu) <= _Tp{0})
	    {
	      gamp = gamma_reciprocal_series(_Tp{1} + mu);
	      gamm = -gamma_reciprocal_series(-mu) / mu;
	    }
	  else
	    {
	      gamp = gamma_reciprocal_series(mu) / mu;
	      gamm = gamma_reciprocal_series(_Tp{1} - mu);
	    }
	  const auto gam1 = (gamm - gamp) / (_Tp{2} * mu);
	  const auto gam2 = (gamm + gamp) / _Tp{2};
	  return gammat_t{mu, gamp, gamm, gam1, gam2};
	}
    }

  /**
   * @brief  Compute the Bessel @f$ J_\nu(x) @f$ and Neumann
   * 	     @f$ N_\nu(x) @f$ functions and their first derivatives
   * 	     @f$ J'_\nu(x) @f$ and @f$ N'_\nu(x) @f$ respectively.
   * 	     These four functions are computed together for numerical
   * 	     stability.
   *
   * @param nu The order of the Bessel functions.
   * @param x  The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    emsr::cyl_bessel_t<_Tp, _Tp, _Tp>
    cyl_bessel_jn_steed(_Tp nu, _Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<_Tp, _Tp, _Tp>;
      const auto s_inf = emsr::infinity(x);
      const auto s_eps = emsr::epsilon(x);
      const auto s_tiny = emsr::lim_min(x);
      const auto s_pi = emsr::pi_v<_Tp>;
      // When the multiplier is N i.e.
      // fp_min = N * min()
      // Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
      //const _Tp s_fp_min = _Tp{20} * emsr::lim_min(nu);
      constexpr int s_max_iter = 15000;
      const auto s_x_min = _Tp{2};
      const auto s_fp_min = emsr::sqrt_min(nu);

      const int n = (x < s_x_min
		    ? std::nearbyint(nu)
		    : std::max(0,
			       static_cast<int>(nu - x + _Tp{1.5L})));

      const auto mu = nu - _Tp(n);
      const auto mu2 = mu * mu;
      const auto xi = _Tp{1} / x;
      const auto xi2 = _Tp{2} * xi;
      const auto _Wronski = xi2 / s_pi;
      int isign = 1;
      auto h = std::max(s_fp_min, nu * xi);
      auto b = xi2 * nu;
      auto d = _Tp{0};
      auto c = h;
      int i;
      for (i = 1; i <= s_max_iter; ++i)
	{
	  b += xi2;
	  d = b - d;
	  if (std::abs(d) < s_fp_min)
	    d = s_fp_min;
	  d = _Tp{1} / d;
	  c = b - _Tp{1} / c;
	  if (std::abs(c) < s_fp_min)
	    c = s_fp_min;
	  const auto del = c * d;
	  h *= del;
	  if (d < _Tp{0})
	    isign = -isign;
	  if (std::abs(del - _Tp{1}) < s_eps)
	    break;
	}
      if (i > s_max_iter)
	return cyl_bessel_jn_asymp(nu, x);

      auto _Jnul = isign * s_fp_min;
      auto _Jpnul = h * _Jnul;
      auto _Jnul1 = _Jnul;
      auto _Jpnu1 = _Jpnul;
      auto fact = nu * xi;
      for (int l = n; l >= 1; --l)
	{
	  const auto _Jnutemp = fact * _Jnul + _Jpnul;
	  fact -= xi;
	  _Jpnul = fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == _Tp{0})
	_Jnul = s_eps;

      const auto f = _Jpnul / _Jnul;
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (x < s_x_min)
	{
	  const auto x2 = x / _Tp{2};
	  const auto pimu = s_pi * mu;
	  const auto fact = (std::abs(pimu) < s_eps
			    ? _Tp{1}
			    : pimu / std::sin(pimu));
	  auto d = -std::log(x2);
	  auto e = mu * d;
	  const auto fact2 = (std::abs(e) < s_eps
			     ? _Tp{1}
			     : std::sinh(e) / e);
	  const auto gamt = gamma_temme(mu);
	  auto ff = (_Tp{2} / s_pi) * fact
		    * (gamt.gamma_1_value * std::cosh(e)
		     + gamt.gamma_2_value * fact2 * d);
	  e = std::exp(e);
	  auto p = e / (s_pi * gamt.gamma_plus_value);
	  auto q = _Tp{1} / (e * s_pi * gamt.gamma_minus_value);
	  const auto pimu2 = pimu / _Tp{2};
	  const auto fact3 = (std::abs(pimu2) < s_eps
			     ? _Tp{1} : std::sin(pimu2) / pimu2 );
	  const auto r = s_pi * pimu2 * fact3 * fact3;
	  auto c = _Tp{1};
	  d = -x2 * x2;
	  auto sum = ff + r * q;
	  auto sum1 = p;
	  int i;
	  for (i = 1; i <= s_max_iter; ++i)
	    {
	      ff = (i * ff + p + q) / (i * i - mu2);
	      c *= d / _Tp(i);
	      p /= _Tp(i) - mu;
	      q /= _Tp(i) + mu;
	      const auto del = c * (ff + r * q);
	      sum += del;
	      const auto del1 = c * p - _Tp(i) * del;
	      sum1 += del1;
	      if (std::abs(del) < s_eps * (_Tp{1} + std::abs(sum)))
		break;
	    }
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_jn_steed: Y-series failed to converge");
	  _Nmu = -sum;
	  _Nnu1 = -sum1 * xi2;
	  _Npmu = mu * xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - f * _Nmu);
	}
      else
	{
	  const auto s_i = std::complex<_Tp>{0, 1};
	  auto a = _Tp{0.25L} - mu2;
	  auto pq = std::complex<_Tp>(-xi / _Tp{2}, _Tp{1});
	  auto b = std::complex<_Tp>(_Tp{2} * x, _Tp{2});
	  auto fact = a * xi / std::norm(pq);
	  auto c = b + s_i * fact * std::conj(pq);
	  auto d = std::conj(b) / std::norm(b);
	  auto dl = c * d;
	  pq *= dl;
	  int i;
	  for (i = 2; i <= s_max_iter; ++i)
	    {
	      a += _Tp(2 * (i - 1));
	      b += s_i * _Tp{2};
	      d = a * d + b;
	      if (std::abs(d) < s_fp_min)
		d = s_fp_min;
	      fact = a / std::norm(c);
	      c = b + fact * std::conj(c);
	      if (std::abs(c) < s_fp_min)
		c = s_fp_min;
	      d = std::conj(d) / std::norm(d);
	      dl = c * d;
	      pq *= dl;
	      if (std::abs(dl - _Tp{1}) < s_eps)
		break;
	    }
	  if (i > s_max_iter)
	    throw std::runtime_error("cyl_bessel_jn_steed: Lentz's method failed");
	  //const auto [p, q] = pq; // This should be a thing.
	  const auto [p, q] = reinterpret_cast<_Tp(&)[2]>(pq);
	  const auto gam = (p - f) / q;
	  _Jmu = std::sqrt(_Wronski / ((p - f) * gam + q));
	  _Jmu = std::copysign(_Jmu, _Jnul);
	  _Nmu = gam * _Jmu;
	  _Npmu = (p + q / gam) * _Nmu;
	  _Nnu1 = mu * xi * _Nmu - _Npmu;
        }
      fact = _Jmu / _Jnul;
      const auto _Jnu = fact * _Jnul1;
      const auto _Jpnu = fact * _Jpnu1;
      if (std::abs(s_pi * x * _Jnu / _Tp{2}) > s_tiny)
	{
	  for (int i = 1; i <= n; ++i)
	    _Nmu = std::exchange(_Nnu1, (mu + i) * xi2 * _Nnu1 - _Nmu);
	  const auto _Nnu = _Nmu;
	  const auto _Npnu = nu * xi * _Nmu - _Nnu1;
	  return bess_t{nu, x, _Jnu, _Jpnu, _Nnu, _Npnu};
	}
      else
	return bess_t{nu, x, _Jnu, _Jpnu, -s_inf, s_inf};
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   * of order @f$ \nu @f$ by various means.
   */
  template<typename _Tp>
    emsr::cyl_bessel_t<_Tp, _Tp, _Tp>
    cyl_bessel_jn(_Tp nu, _Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<_Tp, _Tp, _Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_inf = emsr::infinity(x);
      if (nu < _Tp{0})
	{
	  const auto _Bess = cyl_bessel_jn(-nu, x);
	  const auto sinnupi = emsr::sin_pi(-nu);
	  const auto cosnupi = emsr::cos_pi(-nu);
	  if (std::abs(sinnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(_Tp{1}, cosnupi);
	      return bess_t{nu, x,
			sign * _Bess.J_value, sign * _Bess.J_deriv,
			sign * _Bess.N_value, sign * _Bess.N_deriv};
	    }
	  else if (std::abs(cosnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(_Tp{1}, sinnupi);
	      return bess_t{nu, x,
			-sign * _Bess.N_value, -sign * _Bess.N_deriv,
			 sign * _Bess.J_value,  sign * _Bess.J_deriv};
	    }
	  else
	    {
	      return bess_t{nu, x,
		cosnupi * _Bess.J_value - sinnupi * _Bess.N_value,
		cosnupi * _Bess.J_deriv - sinnupi * _Bess.N_deriv,
		sinnupi * _Bess.J_value + cosnupi * _Bess.N_value,
		sinnupi * _Bess.J_deriv + cosnupi * _Bess.N_deriv};
	    }
	}
      else if (x == _Tp{0})
	{
	  _Tp _Jnu, _Jpnu;
	  if (nu == _Tp{0})
	    {
	      _Jnu = _Tp{1};
	      _Jpnu = _Tp{0};
	    }
	  else if (nu == _Tp{1})
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0.5L};
	    }
	  else
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0};
	    }
	  return bess_t{nu, x, _Jnu, _Jpnu, -s_inf, s_inf};
	}
      else if (x > _Tp{1000})
	return cyl_bessel_jn_asymp(nu, x);
      else
	return cyl_bessel_jn_steed(nu, x);
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   *         of real order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    emsr::cyl_bessel_t<_Tp, _Tp, std::complex<_Tp>>
    cyl_bessel_jn_neg_arg(_Tp nu, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      using bess_t = emsr::cyl_bessel_t<_Tp, _Tp, _Cmplx>;
      constexpr _Cmplx s_i{0, 1};
      if (x >= _Tp{0})
	throw std::domain_error("cyl_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto _Bess = cyl_bessel_jn(nu, -x);
	  const auto phm = emsr::polar_pi(_Tp{1}, -nu);
	  const auto php = emsr::polar_pi(_Tp{1}, nu);
	  const auto cosp = emsr::cos_pi(nu);
	  return bess_t{nu, x,
			  php * _Bess.J_value,
			  -php * _Bess.J_deriv,
			  phm * _Bess.N_value
				+ s_i * _Tp{2} * cosp * _Bess.J_value,
			  -phm * _Bess.N_deriv
				- s_i * _Tp{2} * cosp * _Bess.J_deriv};
	}
    }


  /**
   * @brief  Return the Bessel function of order @f$ \nu @f$:
   * 	     @f$ J_{\nu}(x) @f$.
   *
   * The cylindrical Bessel function is:
   * @f[
   *  J_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  nu  The order of the Bessel function.
   * @param  x   The argument of the Bessel function.
   * @return  The output Bessel function.
   */
  template<typename _Tp>
    _Tp
    cyl_bessel_j(_Tp nu, _Tp x)
    {
      if (x < _Tp{0})
	throw std::domain_error("cyl_bessel_j: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (nu >= _Tp{0} && x * x < _Tp{10} * (nu + _Tp{1}))
	return cyl_bessel_ij_series(nu, x, -1, 200);
      else
	return cyl_bessel_jn(nu, x).J_value;
    }


  /**
   * @brief  Return the Neumann function of order @f$ \nu @f$:
   * 	     @f$ N_{\nu}(x) @f$.
   *
   * The Neumann function is defined by:
   * @f[
   * 	N_{\nu}(x) = \frac{J_{\nu}(x) \cos \nu\pi - J_{-\nu}(x)}
   * 			  {\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   *
   * @param  nu  The order of the Neumann function.
   * @param  x   The argument of the Neumann function.
   * @return  The output Neumann function.
   */
  template<typename _Tp>
    _Tp
    cyl_neumann_n(_Tp nu, _Tp x)
    {
      if (x < _Tp{0})
	throw std::domain_error("cyl_neumann_n: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return cyl_bessel_jn(nu, x).N_value;
    }

  /**
   * @brief  Return the cylindrical Hankel functions of the first
   *         and second kinds and their derivatives.
   *
   */
  template<typename _Tp>
    emsr::cyl_hankel_t<_Tp, _Tp, std::complex<_Tp>>
    cyl_hankel_h1h2(_Tp nu, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx s_i{0, 1};

      _Cmplx ph1 = _Tp{1}, ph2 = _Tp{1};
      if (nu < _Tp{0})
	{
	  ph1 = emsr::polar_pi(_Tp{1}, -nu);
	  ph2 = emsr::polar_pi(_Tp{1}, +nu);
	  nu = -nu;
	}

      // The two _Bess types are different.
      // We might still be able to assign the real output to the complex one.
      if (x < _Tp{0})
	{
	  const auto _Bess = cyl_bessel_jn_neg_arg(nu, x);
	  const auto _H1 = ph1 * (_Bess.J_value + s_i * _Bess.N_value);
	  const auto _dH1 = ph1 * (_Bess.J_deriv + s_i * _Bess.N_deriv);
	  const auto _H2 = ph2 * (_Bess.J_value - s_i * _Bess.N_value);
	  const auto _dH2 = ph2 * (_Bess.J_deriv - s_i * _Bess.N_deriv);
	  return {nu, x, _H1, _dH1, _H2, _dH2};
	}
      else
	{
	  const auto _Bess = cyl_bessel_jn(nu, x);
	  const auto _H1 = ph1 * _Cmplx{_Bess.J_value, _Bess.N_value};
	  const auto _dH1 = ph1 * _Cmplx{_Bess.J_deriv, _Bess.N_deriv};
	  const auto _H2 = ph2 * _Cmplx{_Bess.J_value, -_Bess.N_value};
	  const auto _dH2 = ph2 * _Cmplx{_Bess.J_deriv, -_Bess.N_deriv};
	  return {nu, x, _H1, _dH1, _H2, _dH2};
	}
    }

  /**
   * @brief  Return the cylindrical Hankel function of the first kind
   * 	     @f$ H^{(1)}_\nu(x) @f$.
   *
   * The cylindrical Hankel function of the first kind is defined by:
   * @f[
   *   H^{(1)}_\nu(x) = J_\nu(x) + i N_\nu(x)
   * @f]
   *
   * @param  nu  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    cyl_hankel_1(_Tp nu, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr _Cmplx s_i{0, 1};
      if (nu < _Tp{0})
	return emsr::polar_pi(_Tp{1}, -nu)
	     * cyl_hankel_1(-nu, x);
      else if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < _Tp{0})
	{
	  const auto _Bess = cyl_bessel_jn_neg_arg(nu, x);
	  return _Bess.J_value + s_i * _Bess.N_value;
	}
      else
	{
	  const auto _Bess = cyl_bessel_jn(nu, x);
	  return _Cmplx{_Bess.J_value, _Bess.N_value};
	}
    }


  /**
   *   @brief  Return the cylindrical Hankel function of the second kind
   *           @f$ H^{(2)}_nu(x) @f$.
   *
   *   The cylindrical Hankel function of the second kind is defined by:
   *   @f[
   *     H^{(2)}_\nu(x) = J_\nu(x) - i N_\nu(x)
   *   @f]
   *
   *   @param  nu  The order of the spherical Neumann function.
   *   @param  x  The argument of the spherical Neumann function.
   *   @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    cyl_hankel_2(_Tp nu, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr _Cmplx s_i{0, 1};
      if (nu < _Tp{0})
	return emsr::polar_pi(_Tp{1}, nu)
	     * cyl_hankel_2(-nu, x);
      else if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < _Tp{0})
	{
	  const auto _Bess = cyl_bessel_jn_neg_arg(nu, x);
	  return _Bess.J_value - s_i * _Bess.N_value;
	}
      else
	{
	  const auto _Bess = cyl_bessel_jn(nu, x);
	  return _Cmplx{_Bess.J_value, -_Bess.N_value};
	}
    }


  /**
   * @brief  Compute the spherical Bessel @f$ j_n(x) @f$
   * 	     and Neumann @f$ n_n(x) @f$ functions and their first
   * 	     derivatives @f$ j_n(x) @f$ and @f$ n'_n(x) @f$
   * 	     respectively.
   *
   * @param  n  The order of the spherical Bessel function.
   * @param  x  The argument of the spherical Bessel function.
   * @return  The output derivative of the spherical Neumann function.
   */
  template<typename _Tp>
    emsr::sph_bessel_t<unsigned int, _Tp, _Tp>
    sph_bessel_jn(unsigned int n, _Tp x)
    {
      using bess_t = emsr::sph_bessel_t<unsigned int, _Tp, _Tp>;
      const auto nu = _Tp(n + 0.5L);

      const auto _Bess = cyl_bessel_jn(nu, x);

      const auto factor = (emsr::sqrtpi_v<_Tp> / emsr::sqrt2_v<_Tp>)
			  / std::sqrt(x);

      const auto j_n = factor * _Bess.J_value;
      const auto jp_n = factor * _Bess.J_deriv - j_n / (_Tp{2} * x);
      const auto n_n = factor * _Bess.N_value;
      const auto np_n = factor * _Bess.N_deriv - n_n / (_Tp{2} * x);

      return bess_t{n, x, j_n, jp_n, n_n, np_n};
    }

  /**
   * Return the spherical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    emsr::sph_bessel_t<unsigned int, _Tp, std::complex<_Tp>>
    sph_bessel_jn_neg_arg(unsigned int n, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      using bess_t = emsr::sph_bessel_t<unsigned int, _Tp, _Cmplx>;
      if (x >= _Tp{0})
	throw std::domain_error("sph_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto nu = _Tp(n + 0.5L);
	  const auto _Bess = cyl_bessel_jn_neg_arg(nu, x);

	  const auto factor
	    = (emsr::sqrtpi_v<_Tp> / emsr::sqrt2_v<_Tp>)
	      / std::sqrt(_Cmplx(x));

	  const auto j_n = factor * _Bess.J_value;
	  const auto jp_n = factor * _Bess.J_deriv
			    - j_n / (_Tp{2} * x);
	  const auto n_n = factor * _Bess.N_value;
	  const auto np_n = factor * _Bess.N_deriv
			    - n_n / (_Tp{2} * x);

	  return bess_t{n, x, j_n, jp_n, n_n, np_n};
	}
    }


  /**
   * @brief  Return the spherical Bessel function @f$ j_n(x) @f$ of order n
   * and non-negative real argument @c x.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *   j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @param  n  The non-negative integral order
   * @param  x  The non-negative real argument
   * @return  The output spherical Bessel function.
   */
  template<typename _Tp>
    _Tp
    sph_bessel(unsigned int n, _Tp x)
    {
      if (x < _Tp{0})
	throw std::domain_error("sph_bessel: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == _Tp{0})
	{
	  if (n == 0)
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	return sph_bessel_jn(n, x).j_value;
    }


  /**
   * @brief  Return the spherical Neumann function @f$ n_n(x) @f$ of order n
   * and non-negative real argument @c x.
   *
   * The spherical Neumann function is defined by:
   * @f[
   *  n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @param  n  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    _Tp
    sph_neumann(unsigned int n, _Tp x)
    {
      if (x < _Tp{0})
	throw std::domain_error("sph_neumann: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == _Tp{0})
	return -emsr::infinity(x);
      else
	return sph_bessel_jn(n, x).n_value;
    }


  /**
   * @brief  Return the spherical Hankel function of the first kind
   * 	     @f$ h^{(1)}_n(x) @f$.
   *
   * The spherical Hankel function of the first kind is defined by:
   * @f[
   *   h^{(1)}_n(x) = j_n(x) + i n_n(x)
   * @f]
   *
   * @param  n  The order of the spherical Neumann function.
   * @param  x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    sph_hankel_1(unsigned int n, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < _Tp{0})
	{
	  const auto _Bess = sph_bessel_jn_neg_arg(n, x);
	  return _Bess.j_value + s_i * _Bess.n_value;
	}
      else
	{
	  const auto _Bess = sph_bessel_jn(n, x);
	  return _Cmplx{_Bess.j_value, _Bess.n_value};
	}
    }


  /**
   * @brief  Return the spherical Hankel function of the second kind
   * 	     @f$ h^{(2)}_n(x) @f$.
   *
   * The spherical Hankel function of the second kind is defined by:
   * @f[
   *   h^{(2)}_n(x) = j_n(x) - i n_n(x)
   * @f]
   *
   * @param  n  The non-negative integral order
   * @param  x  The non-negative real argument
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    sph_hankel_2(unsigned int n, _Tp x)
    {
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < _Tp{0})
	{
	  const auto _Bess = sph_bessel_jn_neg_arg(n, x);
	  return _Bess.j_value - s_i * _Bess.n_value;
	}
      else
	{
	  const auto _Bess = sph_bessel_jn(n, x);
	  return _Cmplx{_Bess.j_value, -_Bess.n_value};
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_BESSEL_TCC
