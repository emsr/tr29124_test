// Special functions -*- C++ -*-

// Copyright (C) 2006-2019 Free Software Foundation, Inc.
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

#include <emsr/continued_fractions.h>

#define BESSEL_DEBUG 1
#ifdef BESSEL_DEBUG
#  include <iostream>
#  include <iomanip>
#endif
#include <utility> // For exchange.

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
  template<typename _Tnu, typename Tp>
    constexpr Tp
    cyl_bessel_ij_series(_Tnu nu, Tp x, int sgn,
			   unsigned int max_iter)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon<_Real>();
      if (std::abs(x) < s_eps)
	{
	  if (nu == _Tnu{0})
	    return Tp{1};
	  else
	    return Tp{0};
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
  template<typename _Tnu, typename Tp>
    struct cyl_bessel_asymp_sums_t
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
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
  template<typename _Tnu, typename Tp>
    constexpr cyl_bessel_asymp_sums_t<_Tnu, Tp>
    cyl_bessel_asymp_sums(_Tnu nu, Tp x, int sgn)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using bess_t = cyl_bessel_asymp_sums_t<_Tnu, Tp>;
      const auto s_eps = emsr::epsilon<_Real>();
      const auto __2nu = _Real{2} * nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto r8x = Tp{1} / (_Real{8} * x);
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
  template<typename _Tnu, typename Tp>
    constexpr emsr::cyl_bessel_t<_Tnu, Tp, Tp>
    cyl_bessel_jn_asymp(_Tnu nu, Tp x)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using bess_t = emsr::cyl_bessel_t<_Tnu, Tp, Tp>;
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
   * Compute ratios of Bessel functions using the S-fraction.
   *
   * @param  nu  The order of the Hankel ratio.
   * @param  x   The argument of the Hankel ratio.
   * @param  zeta A variable encapsulating the regular and irregular
   *           Bessel functions; Use @f$ \zeta = (iz)^2 @f$ for @f$ J_\nu(z) @f$
   *           and @f$ \zeta = z^2 @f$ for @f$ I_\nu(z) @f$.
   */
  template<typename _Tnu, typename Tp, typename _Tzeta>
    std::complex<emsr::num_traits_t<
		 emsr::fp_promote_t<_Tnu, Tp, _Tzeta>>>
    cyl_bessel_ratio_s_frac(_Tnu nu, Tp z, _Tzeta zeta)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      auto a_J
	= [nu, z, zeta](std::size_t k, Tp)
	  {
	    using type = decltype(_Tnu{} * _Tzeta{});
	    if (k == 1)
	      return type(z / (_Tnu{2} * nu + _Tnu{2}));
	    else
	      return zeta
		   / (_Tnu{4} * (nu + _Tnu(k - 1)) * (nu + _Tnu(k)));
	  };
      using _AFun = decltype(a_J);

      auto b_J = [](std::size_t, Tp) -> _Real { return _Real{1}; };
      using _BFun = decltype(b_J);

      auto w_J = [](std::size_t, Tp) -> _Real { return _Real{0}; };
      using _WFun = decltype(w_J);

      emsr::SteedContinuedFraction<Tp, _AFun, _BFun, _WFun>
      _J(a_J, b_J, w_J);

      // b_0 is 0 not 1 so subtract 1.
      return _J(z) - _Real{1};
    }

  /**
   * 
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_j_ratio_s_frac(_Tnu nu, Tp z)
    {
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto iz = _Cmplx{0, 1} * z;
      const auto zeta = iz * iz;
      const auto _Jrat = cyl_bessel_ratio_s_frac(nu, z, zeta);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Jrat);
      else
	return _Jrat;
    }

  /**
   * 
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_i_ratio_s_frac(_Tnu nu, Tp z)
    {
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto zeta = z * z;
      const auto _Irat = cyl_bessel_ratio_s_frac(nu, z, zeta);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Irat);
      else
	return _Irat;
    }

  /**
   * Compute ratios of Hankel functions using the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    std::complex<emsr::num_traits_t<
		 emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_ratio_j_frac(_Tnu nu, Tp z, Tp sgn)
    {
      using _Val = emsr::fp_promote_t<_Tnu, Tp>;
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto zeta = _Cmplx{0, 2} * z;
      using _Tzeta = decltype(zeta);

      auto a_H
	= [nu](std::size_t k, Tp)
	  {
	    const auto kk = _Tnu(2 * k - 1) / _Tnu{2};
	    return (nu - kk) * (nu + kk);
	  };
      using _NumFun = decltype(a_H);

      auto b_H
	= [zeta, sgn](std::size_t k, Tp)
	  { return sgn * _Tzeta(2 * k) + zeta; };
      using _DenFun = decltype(b_H);

      auto w_H
	= [zeta](std::size_t k, Tp)
	  { return _Tzeta(k) + zeta / _Tzeta{2}; };
      using _TailFun = decltype(w_H);

      emsr::SteedContinuedFraction<Tp, _NumFun, _DenFun, _TailFun>
      _H(a_H, b_H, w_H);

      return (_Tzeta(2 * nu + 1) + sgn * zeta) / (Tp{2} * z)
	   + sgn * (_H(z) - b_H(0, Tp{})) / z;
    }

  /**
   * Return the Hankel function ratio of the first kind from the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    inline std::complex<emsr::num_traits_t<
			emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_1_ratio_j_frac(_Tnu nu, Tp z)
    { return cyl_hankel_ratio_j_frac(nu, z, Tp{-1}); }

  /**
   * Return the Hankel function ratio of the second kind from the J-fraction.
   */
  template<typename _Tnu, typename Tp>
    inline std::complex<emsr::num_traits_t<
			emsr::fp_promote_t<_Tnu, Tp>>>
    cyl_hankel_2_ratio_j_frac(_Tnu nu, Tp z)
    { return cyl_hankel_ratio_j_frac(nu, z, Tp{+1}); }

  /**
   * Return the modified Bessel function ratio of the second kind
   * from the J-fraction ratios of Hankel functions.
   */
  template<typename _Tnu, typename Tp,
	   typename _Val = emsr::fp_promote_t<_Tnu, Tp>>
    std::conditional_t<emsr::is_complex_v<_Val>,
			std::complex<emsr::num_traits_t<_Val>>,
			_Val>
    cyl_bessel_k_ratio_j_frac(_Tnu nu, Tp z)
    {
      using _Real = emsr::num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto s_i = _Cmplx{0, 1};
      const auto s_pi = emsr::pi_v<_Real>;
      const auto ph = std::arg(z);

      _Cmplx _Krat;
      if (ph > -s_pi && ph <= s_pi / _Real{2})
	_Krat = s_i * cyl_hankel_1_ratio_j_frac(nu, s_i * z);
      else
	_Krat = -s_i * cyl_hankel_2_ratio_j_frac(nu, -s_i * z);

      if constexpr (!emsr::is_complex_v<_Val>)
	return std::real(_Krat);
      else
	return _Krat;
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
  template<typename Tp>
    emsr::gamma_temme_t<Tp>
    gamma_temme(Tp mu)
    {
      using gammat_t = emsr::gamma_temme_t<Tp>;
      const auto s_eps = emsr::epsilon(mu);
      const auto s_gamma_E = emsr::egamma_v<Tp>;

      if (std::abs(mu) < s_eps)
	return gammat_t{mu, Tp{1}, Tp{1}, -s_gamma_E, Tp{1}};
      else
	{
	  Tp gamp, gamm;
	  if (std::real(mu) <= Tp{0})
	    {
	      gamp = gamma_reciprocal_series(Tp{1} + mu);
	      gamm = -gamma_reciprocal_series(-mu) / mu;
	    }
	  else
	    {
	      gamp = gamma_reciprocal_series(mu) / mu;
	      gamm = gamma_reciprocal_series(Tp{1} - mu);
	    }
	  const auto gam1 = (gamm - gamp) / (Tp{2} * mu);
	  const auto gam2 = (gamm + gamp) / Tp{2};
	  return gammat_t{mu, gamp, gamm, gam1, gam2};
	}
    }

  /**
   * A little return type for the Temme series.
   */
  template<typename Tp>
    struct bessel_nk_series_t
    {
      Tp _Z_mu;
      Tp _Z_mup1;
    };

  /**
   * This routine computes the dominant cylindrical bessel function solutions
   * by series summation for order @f$ |\mu| < 1/2 @f$.
   *
   * @param mu The order of the Bessel functions @f$ |\mu| < 1/2 @f$.
   * @param x  The argument of the Bessel functions.
   * @param modified If true solve for the modified Bessel function
   *                   @f$ K_\mu @f$, otherwise solve for the Neumann function
   *                   @f$ N_\mu @f$,
   * @return A structure containing Z_{\mu} and Z_{\mu+1}.
   */
  template<typename Tp>
    bessel_nk_series_t<Tp>
    cyl_bessel_nk_series(Tp mu, Tp x, bool modified = false,
			   int max_iter = 100)
    {
      const auto s_eps = emsr::epsilon<Tp>();
      const auto s_pi = emsr::pi_v<Tp>;
      const auto xi = Tp{1} / x;
      const auto x2 = x / Tp{2};

      const auto fact = Tp{1} / emsr::detail::sinc_pi(mu);
      const auto lx2 = -std::log(x2);
      const auto arg = mu * lx2;
      const auto fact2 = emsr::detail::sinhc(arg);
      const auto gamt = emsr::detail::gamma_temme(mu);
      const auto norm = modified ? Tp{-1} : Tp{2} / s_pi;
      auto ff = norm * fact
		* (gamt.gamma_1_value * std::cosh(arg)
		 + gamt.gamma_2_value * fact2 * lx2);
      const auto e = std::exp(arg);
      auto p = norm * e / (Tp{2} * gamt.gamma_plus_value);
      auto q = norm / (e * Tp{2} * gamt.gamma_minus_value);
      const auto fact3 = modified
			 ? Tp{0}
			 : emsr::detail::sinc_pi(mu / Tp{2});
      const auto r = modified
		     ? Tp{0}
		     : fact3 * fact3 * s_pi * s_pi * mu / Tp{2};
      auto c = Tp{1};
      const auto d = modified ? x2 * x2 : -x2 * x2;
      auto sum_mu = ff + r * q;
      auto sum_mup1 = p;
      int i;
      for (i = 1; i <= max_iter; ++i)
	{
	  ff = (i * ff + p + q)
	       / ((Tp(i) - mu) * (Tp(i) + mu));
	  c *= d / Tp(i);
	  p /= Tp(i) - mu;
	  q /= Tp(i) + mu;
	  const auto del_mu = c * (ff + r * q);
	  sum_mu += del_mu;
	  const auto del_mup1 = c * p - Tp(i) * del_mu;
	  sum_mup1 += del_mup1;
	  if (std::abs(del_mu) < s_eps * std::abs(sum_mu))
	    break;
	}
      if (i > max_iter)
	throw std::runtime_error("cyl_bessel_nk_series: Series failed to converge");
      auto _N_mu = -sum_mu;
      auto _N_mup1 = -Tp{2} * xi * sum_mup1;

      return {_N_mu, _N_mup1};
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
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, Tp>
    cyl_bessel_jn_steed(Tp nu, Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Tp>;
      const auto s_inf = emsr::infinity(x);
      const auto s_eps = emsr::epsilon(x);
      const auto s_tiny = emsr::lim_min(x);
      const auto s_pi = emsr::pi_v<Tp>;
      // When the multiplier is N i.e.
      // fp_min = N * min()
      // Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
      //const Tp s_fp_min = Tp{20} * emsr::lim_min(nu);
      constexpr int s_max_iter = 15000;
      const auto s_x_min = Tp{2};
      const auto s_fp_min = emsr::sqrt_min(nu);

      const int n = (x < s_x_min
		    ? std::nearbyint(nu)
		    : std::max(0,
			       static_cast<int>(nu - x + Tp{1.5L})));

      const auto xi = Tp{1} / x;
      const auto xi2 = Tp{2} * xi;
      const auto _Wronski = xi2 / s_pi;
      int isign = 1;
//#if 1
      auto h = std::max(s_fp_min, nu * xi);
      auto b = xi2 * nu;
      auto d = Tp{0};
      auto c = h;
      int i;
      for (i = 1; i <= s_max_iter; ++i)
	{
	  b += xi2;
	  d = b - d;
	  if (std::abs(d) < s_fp_min)
	    d = s_fp_min;
	  d = Tp{1} / d;
	  c = b - Tp{1} / c;
	  if (std::abs(c) < s_fp_min)
	    c = s_fp_min;
	  const auto del = c * d;
	  h *= del;
	  if (d < Tp{0})
	    isign = -isign;
	  if (std::abs(del - Tp{1}) < s_eps)
	    break;
	}
      if (i > s_max_iter)
	return cyl_bessel_jn_asymp(nu, x);
//#else
      Tp h2;
      try
	{
	  h2 = nu * xi - cyl_bessel_j_ratio_s_frac(nu, x);
	}
      catch (...)
	{
	  return cyl_bessel_jn_asymp(nu, x);
	}
//#endif
#ifdef BESSEL_DEBUG
std::cerr << ' ' << std::setw(10) << nu
	  << ' ' << std::setw(10) << x
	  << ' ' << std::setw(2) << isign
	  << ' ' << std::setw(10) << h2
	  << ' ' << std::setw(10) << h
	  << ' ' << std::setw(14) << h2 - h
	  << '\n';
#endif

      // Recur recessive solution downward to |mu| < 1/2.
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
      if (_Jnul == Tp{0})
	_Jnul = s_eps;

      const auto mu = nu - Tp(n);
      const auto mu2 = mu * mu;
      const auto f = _Jpnul / _Jnul;
      Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (x < s_x_min)
	{
	  const auto _Z = cyl_bessel_nk_series(mu, x);
	  _Nmu = _Z._Z_mu;
	  _Nnu1 = _Z._Z_mup1;
	  _Npmu = mu * xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - f * _Nmu);
	}
      else
	{
	  auto pq = cyl_hankel_1_ratio_j_frac(mu, x);
	  auto [p, q] = reinterpret_cast<Tp(&)[2]>(pq);
	  // FIXME: Is this right?
	  q = std::abs(q);
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
      if (std::abs(s_pi * x * _Jnu / Tp{2}) > s_tiny)
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
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, Tp>
    cyl_bessel_jn(Tp nu, Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_inf = emsr::infinity(x);
      const auto s_pi = emsr::pi_v<Tp>;
      if (nu < Tp{0})
	{
	  const auto _Bess = cyl_bessel_jn(-nu, x);
	  const auto sinnupi = sin_pi(-nu);
	  const auto cosnupi = cos_pi(-nu);
	  if (std::abs(sinnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(Tp{1}, cosnupi);
	      return bess_t{nu, x,
			sign * _Bess.J_value, sign * _Bess.J_deriv,
			sign * _Bess.N_value, sign * _Bess.N_deriv};
	    }
	  else if (std::abs(cosnupi) < s_eps)
	    { // Carefully preserve +-inf.
	      const auto sign = std::copysign(Tp{1}, sinnupi);
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
      else if (x == Tp{0})
	{
	  Tp _Jnu, _Jpnu;
	  if (nu == Tp{0})
	    {
	      _Jnu = Tp{1};
	      _Jpnu = Tp{0};
	    }
	  else if (nu == Tp{1})
	    {
	      _Jnu = Tp{0};
	      _Jpnu = Tp{0.5L};
	    }
	  else
	    {
	      _Jnu = Tp{0};
	      _Jpnu = Tp{0};
	    }
	  return bess_t{nu, x, _Jnu, _Jpnu, -s_inf, s_inf};
	}
      else if (x > Tp{1000})
	return cyl_bessel_jn_asymp(nu, x);
      else
	return cyl_bessel_jn_steed(nu, x);
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   *         of real order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename Tp>
    emsr::cyl_bessel_t<Tp, Tp, std::complex<Tp>>
    cyl_bessel_jn_neg_arg(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, _Cmplx>;
      constexpr _Cmplx s_i{0, 1};
      if (x >= Tp{0})
	throw std::domain_error("cyl_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto _Bess = cyl_bessel_jn(nu, -x);
	  const auto phm = polar_pi(Tp{1}, -nu);
	  const auto php = polar_pi(Tp{1}, nu);
	  const auto cosp = cos_pi(nu);
	  return bess_t{nu, x,
			  php * _Bess.J_value,
			  -php * _Bess.J_deriv,
			  phm * _Bess.N_value
				+ s_i * Tp{2} * cosp * _Bess.J_value,
			  -phm * _Bess.N_deriv
				- s_i * Tp{2} * cosp * _Bess.J_deriv};
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
  template<typename Tp>
    Tp
    cyl_bessel_j(Tp nu, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("cyl_bessel_j: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (nu >= Tp{0} && x * x < Tp{10} * (nu + Tp{1}))
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
  template<typename Tp>
    Tp
    cyl_neumann_n(Tp nu, Tp x)
    {
      if (x < Tp{0})
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
  template<typename Tp>
    emsr::cyl_hankel_t<Tp, Tp, std::complex<Tp>>
    cyl_hankel_h1h2(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      constexpr _Cmplx s_i{0, 1};

      _Cmplx ph1 = Tp{1}, ph2 = Tp{1};
      if (nu < Tp{0})
	{
	  ph1 = polar_pi(Tp{1}, -nu);
	  ph2 = polar_pi(Tp{1}, +nu);
	  nu = -nu;
	}

      // The two _Bess types are different.
      // We might still be able to assign the real output to the complex one.
      if (x < Tp{0})
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
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_1(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr _Cmplx s_i{0, 1};
      if (nu < Tp{0})
	return polar_pi(Tp{1}, -nu)
	     * cyl_hankel_1(-nu, x);
      else if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
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
  template<typename Tp>
    std::complex<Tp>
    cyl_hankel_2(Tp nu, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      const auto s_nan = emsr::quiet_NaN(x);
      constexpr _Cmplx s_i{0, 1};
      if (nu < Tp{0})
	return polar_pi(Tp{1}, nu)
	     * cyl_hankel_2(-nu, x);
      else if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
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
  template<typename Tp>
    emsr::sph_bessel_t<unsigned int, Tp, Tp>
    sph_bessel_jn(unsigned int n, Tp x)
    {
      using bess_t = emsr::sph_bessel_t<unsigned int, Tp, Tp>;
      const auto nu = Tp(n + 0.5L);

      const auto _Bess = cyl_bessel_jn(nu, x);

      const auto factor = (emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>)
			  / std::sqrt(x);

      const auto j_n = factor * _Bess.J_value;
      const auto jp_n = factor * _Bess.J_deriv - j_n / (Tp{2} * x);
      const auto n_n = factor * _Bess.N_value;
      const auto np_n = factor * _Bess.N_deriv - n_n / (Tp{2} * x);

      return bess_t{n, x, j_n, jp_n, n_n, np_n};
    }

  /**
   * Return the spherical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename Tp>
    emsr::sph_bessel_t<unsigned int, Tp, std::complex<Tp>>
    sph_bessel_jn_neg_arg(unsigned int n, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      using bess_t
	= emsr::sph_bessel_t<unsigned int, Tp, _Cmplx>;
      if (x >= Tp{0})
	throw std::domain_error("sph_bessel_jn_neg_arg: non-negative argument");
      else
	{
	  const auto nu = Tp(n + 0.5L);
	  const auto _Bess = cyl_bessel_jn_neg_arg(nu, x);

	  const auto factor
	    = (emsr::sqrtpi_v<Tp> / emsr::sqrt2_v<Tp>)
	      / std::sqrt(_Cmplx(x));

	  const auto j_n = factor * _Bess.J_value;
	  const auto jp_n = factor * _Bess.J_deriv
			    - j_n / (Tp{2} * x);
	  const auto n_n = factor * _Bess.N_value;
	  const auto np_n = factor * _Bess.N_deriv
			    - n_n / (Tp{2} * x);

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
  template<typename Tp>
    Tp
    sph_bessel(unsigned int n, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("sph_bessel: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == Tp{0})
	{
	  if (n == 0)
	    return Tp{1};
	  else
	    return Tp{0};
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
  template<typename Tp>
    Tp
    sph_neumann(unsigned int n, Tp x)
    {
      if (x < Tp{0})
	throw std::domain_error("sph_neumann: bad argument");
      else if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x == Tp{0})
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
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_1(unsigned int n, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      constexpr _Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
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
  template<typename Tp>
    std::complex<Tp>
    sph_hankel_2(unsigned int n, Tp x)
    {
      using _Cmplx = std::complex<Tp>;
      constexpr _Cmplx s_i{0, 1};
      const auto s_nan = emsr::quiet_NaN(x);
      if (std::isnan(x))
	return _Cmplx{s_nan, s_nan};
      else if (x < Tp{0})
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
} // namespace std

#endif // SF_BESSEL_TCC
