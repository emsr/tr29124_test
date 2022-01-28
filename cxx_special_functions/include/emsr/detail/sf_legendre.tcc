
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

/** @file bits/sf_legendre.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1) https://dlmf.nist.gov/14 - Chapter 14 Legendre and Related Functions
// (2) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 8, pp. 331-341
// (3) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (4) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 252-254
// (5) Computation of Special Functions, Shanjie Zhang, Jianming Jin,
//     John Wiley & Sons (1996)

#ifndef SF_LEGENDRE_TCC
#define SF_LEGENDRE_TCC 1

#include <vector>
#include <complex>

#include <emsr/sf_hyperg.h>
#include <emsr/specfun_state.h>
#include <emsr/math_constants.h>
#include <emsr/quadrature_point.h>

namespace emsr
{
namespace detail
{
  /**
   * @brief  Return the Legendre polynomial by upward recursion
   * 	     on degree @f$ l @f$.
   *
   * The Legendre function of degree @f$ l @f$ and argument @f$ x @f$,
   * @f$ P_l(x) @f$, is defined by:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   * This can be expressed as a series:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\sum_{k=0}^{\lfloor l/2 \rfloor}
   *            \frac{(-1)^k(2l-2k)!}{k!(l-k)!(l-2k)!}x^{l-2k}
   * @f]
   *
   * @param  l  The degree of the Legendre polynomial.  @f$ l >= 0 @f$.
   * @param  x  The argument of the Legendre polynomial.
   */
  template<typename Tp>
    emsr::legendre_p_t<Tp>
    legendre_p(unsigned int l, Tp x)
    {
      using _Real = emsr::num_traits_t<Tp>;
      using ret_t = emsr::legendre_p_t<Tp>;

      const auto lge1 = l >= 1 ? Tp{+1} : Tp{0};
      const auto lge2 = l >= 2 ? Tp{+1} : Tp{0};
      const auto s_NaN = emsr::quiet_NaN(_Real{});

      if (std::isnan(x))
	return {l, s_NaN, s_NaN, s_NaN, s_NaN};
      else if (x == _Real{+1})
	return {l, x, Tp{+1}, lge1, lge2};
      else if (x == _Real{-1})
	return l % 2 == 1
		? ret_t{l, x, Tp{-1}, +lge1, -lge2}
		: ret_t{l, x, Tp{+1}, -lge1, +lge2};
      else
	{
	  auto _P_lm2 = Tp{1};
	  if (l == 0)
	    return {l, x, _P_lm2, Tp{0}, Tp{0}};

	  auto _P_lm1 = x;
	  if (l == 1)
	    return {l, x, _P_lm1, _P_lm2, Tp{0}};

	  auto _P_l = Tp{2} * x * _P_lm1 - _P_lm2
		    - (x * _P_lm1 - _P_lm2) / Tp{2};
	  for (unsigned int ll = 3; ll <= l; ++ll)
	    {
	      _P_lm2 = _P_lm1;
	      _P_lm1 = _P_l;
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      _P_l = Tp{2} * x * _P_lm1 - _P_lm2
		    - (x * _P_lm1 - _P_lm2) / Tp(ll);
	    }

	  return {l, x, _P_l, _P_lm1, _P_lm2};
	}
    }

  /**
   * Legendre q series.
   */
  template<typename Tp>
    Tp
    legendre_q_series(unsigned int l, Tp x)
    {
      const auto num1 = Tp(l + 1) / Tp{2};
      const auto num2 = Tp(l + 2) / Tp{2};
      const auto den = Tp(2 * l + 3) / Tp{2};
      const auto rx = Tp{1} / x;
      const auto rx2 = rx * rx;

      const auto sum = emsr::hyperg(num1, num2, den, rx2);

      auto fact = rx;
      for (unsigned int k = 1; k <= l; ++k)
	fact *= Tp(k) * rx / Tp(2 * k - 1);

      return fact * sum;
    }

  /**
   * @brief Return the Legendre function of the second kind
   *        by upward recursion on degree @f$ l @f$.
   *
   * The Legendre function of the second kind of degree @f$ l @f$
   * and argument @f$ x @f$, @f$ Q_l(x) @f$, is defined by:
   * @f[
   *   Q_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   *
   * @param l The degree of the Legendre function.  @f$ l >= 0 @f$.
   * @param x The argument of the Legendre function.  @f$ |x| <= 1 @f$.
   */
  template<typename Tp>
    emsr::legendre_q_t<Tp>
    legendre_q(unsigned int l, Tp x)
    {
      using _Real = emsr::num_traits_t<Tp>;
      const auto s_eps = emsr::epsilon(_Real{});
      const auto s_inf = emsr::infinity(_Real{});
      if (std::isnan(x))
	{
	  const auto s_NaN = emsr::quiet_NaN(_Real{});
	  return {l, x, s_NaN, s_NaN, s_NaN};
	}
      else if (std::abs(x - _Real{1}) < s_eps)
	return {l, x, s_inf, s_inf, s_inf};
      else if (std::abs(x + _Real{1}) < s_eps)
	{
	  const auto sgn = (l & 1 ? +1 : -1);
	  return {l, x, sgn * s_inf, -sgn * s_inf, sgn * s_inf};
	}
      else if (std::abs(x) < _Real{1})
	{
	  auto _Q_lm2 = Tp{0.5L} * std::log((Tp{1} + x) / (Tp{1} - x));
	  if (l == 0)
	    return {l, x, _Q_lm2, Tp{0}, Tp{0}};
	  auto _Q_lm1 = x * _Q_lm2 - Tp{1};
	  if (l == 1)
	    return {l, x, _Q_lm1, _Q_lm2, Tp{0}};
	  auto _Q_l = Tp{2} * x * _Q_lm1 - _Q_lm2
		    - (x * _Q_lm1 - _Q_lm2) / Tp{2};
	  for (unsigned int ll = 3; ll <= l; ++ll)
	    {
	      _Q_lm2 = _Q_lm1;
	      _Q_lm1 = _Q_l;
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      // Is this still true for Q?
	      _Q_l = Tp{2} * x * _Q_lm1 - _Q_lm2
		    - (x * _Q_lm1 - _Q_lm2) / Tp(ll);
	    }

	  return {l, x, _Q_l, _Q_lm1, _Q_lm2};
	}
      else
	{
          // FIXME: Knock these out with one series.
	  const auto _Q_l = legendre_q_series(l, x);
	  if (l == 0)
	    return {l, x, _Q_l, Tp{0}, Tp{0}};
	  const auto _Q_lm1 = legendre_q_series(l - 1, x);
	  if (l == 1)
	    return {l, x, _Q_l, _Q_lm1, Tp{0}};
	  const auto _Q_lm2 = (Tp(2 * l - 1) * x * _Q_lm1
			     - Tp(l) * _Q_l) / Tp(l - 1);
	  return {l, x, _Q_l, _Q_lm1, _Q_lm2};
	}
    }

  /**
   * @brief  Return the associated Legendre function by recursion
   * 	     on @f$ l @f$ and downward recursion on m.
   *
   * The associated Legendre function is derived from the Legendre function
   * @f$ P_l(x) @f$ by the Rodrigues formula:
   * @f[
   *   P_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
   * @f]
   * @note The Condon-Shortley phase factor @f$ (-1)^m @f$ is absent
   *       by default.
   * @note @f$ P_l^m(x) = 0 @f$ if @f$ m > l @f$.
   *
   * @param  l  The degree of the associated Legendre function.
   * 		@f$ l >= 0 @f$.
   * @param  m  The order of the associated Legendre function.
   * @param  x  The argument of the associated Legendre function.
   * @param  phase  The phase of the associated Legendre function.
   *                  Use -1 for the Condon-Shortley phase convention.
   */
  template<typename Tp>
    emsr::assoc_legendre_p_t<Tp>
    assoc_legendre_p(unsigned int l, unsigned int m, Tp x,
		       Tp phase = Tp{+1})
    {
      using _Real = emsr::num_traits_t<Tp>;
      if (m > l)
	return {l, m, x, Tp{0}, Tp{0}, Tp{0}};
      else if (std::isnan(x))
	{
	  const auto _NaN = emsr::quiet_NaN(_Real{});
	  return {l, m, x, _NaN, _NaN, _NaN, phase};
	}
      else if (m == 0)
	{
	  const auto _P_l = legendre_p(l, x);
	  return {l, m, x, _P_l.P_l, _P_l.P_lm1, _P_l.P_lm2, phase};
	}
      else
	{
	  auto _P_mm = Tp{1};
	  if (m > 0)
	    {
	      // Two square roots seem more accurate more of the time
	      // than just one.
	      auto root = std::sqrt(Tp{1} - x) * std::sqrt(Tp{1} + x);
	      auto fact = Tp{1};
	      for (unsigned int i = 1; i <= m; ++i)
		{
		  // N. B. Condon-Shortley would use phase = -1 here.
		  _P_mm *= phase * fact * root;
		  fact += Tp{2};
		}
	    }
	  if (l == m)
	    return {l, m, x, _P_mm, Tp{0}, Tp{0}, phase};

	  auto _P_mp1m = Tp(2 * m + 1) * x * _P_mm;
	  if (l == m + 1)
	    return {l, m, x, _P_mp1m, _P_mm, Tp{0}, phase};

	  auto _P_lm2m = _P_mm;
	  auto _P_lm1m = _P_mp1m;
	  auto _P_lm = (Tp(2 * m + 3) * x * _P_lm1m
		     - Tp(2 * m + 1) * _P_lm2m) / Tp{2};
	  for (unsigned int j = m + 3; j <= l; ++j)
	    {
	      _P_lm2m = _P_lm1m;
	      _P_lm1m = _P_lm;
	      _P_lm = (Tp(2 * j - 1) * x * _P_lm1m
		      - Tp(j + m - 1) * _P_lm2m) / Tp(j - m);
	    }

	  return {l, m, x, _P_lm, _P_lm1m, _P_lm2m};
	}
    }

  template<typename Tp>
    emsr::assoc_legendre_q_t<Tp>
    assoc_legendre_q(unsigned int l, unsigned int m, Tp x,
		       Tp phase = Tp{+1})
    {
      using _Real = emsr::num_traits_t<Tp>;
      if (std::isnan(x))
	{
	  const auto _NaN = emsr::quiet_NaN(_Real{});
	  return {l, m, x, _NaN, _NaN, _NaN, phase};
	}
      else if (std::abs(x) < _Real{1})
	{
	  // Find Q_l^0 and Q_l^1 by upward recurrence on l.
	  const auto fact = (Tp{1} - x) * (Tp{1} + x);
	  const auto root = std::sqrt(Tp{1} - x) * std::sqrt(Tp{1} + x);
	  const auto log = std::log((Tp{1} + x) / (Tp{1} - x)) / Tp{2};

	  const auto _Q_00 = log;
	  const auto _Q_01 = phase / root;
	  if (l == 0)
	    {
	      if (m == 0)
		return {l, m, x, _Q_00, Tp{0}, Tp{0}, phase}; // FIXME?
	      else if (m == 1)
		return {l, m, x, _Q_01, _Q_00, Tp{0}, phase}; // FIXME?
	    }

	  const auto _Q_10 = x * _Q_00 - Tp{1};
	  const auto _Q_11 = phase * root * (log + x / fact);
	  if (l == 1)
	    {
	      if (m == 0)
		return {l, m, x, _Q_10, Tp{0}, Tp{0}, phase}; // FIXME?
	      else if (m == 1)
		return {l, m, x, _Q_11, _Q_10, Tp{0}, phase}; // FIXME?
	    }

	  auto _Q_lm20 = _Q_00;
	  auto _Q_lm10 = _Q_10;
	  auto _Q_lm21 = _Q_01;
	  auto _Q_lm11 = _Q_11;
	  auto _Q_l0 = (Tp{3} * x * _Q_lm10 - Tp{1} * _Q_lm20) / Tp{2};
	  auto _Q_l1 = (Tp{3} * x * _Q_lm11 - Tp{2} * _Q_lm21);
	  for (unsigned int k = 3; k <= l; ++k)
	    {
	      _Q_lm20 = _Q_lm10;
	      _Q_lm10 = _Q_l0;
	      _Q_l0 = (Tp(2 * k - 1) * x * _Q_lm10
		    - Tp(k - 1) * _Q_lm20) / Tp(k);

	      _Q_lm21 = _Q_lm11;
	      _Q_lm11 = _Q_l1;
	      _Q_l1 = (Tp(2 * k - 1) * x * _Q_lm11
		    - Tp(k) * _Q_lm21) / Tp(k - 1);
	    }
	  if (m == 0)
	    return {l, m, x, _Q_l0, Tp{0}, Tp{0}, phase}; // FIXME?
	  else if (m == 1)
	    return {l, m, x, _Q_l1, _Q_l0, Tp{0}, phase}; // FIXME?

	  // Find Q_l^m by upward recurrence on m.
	  auto _Q_lmm2 = _Q_l0;
	  auto _Q_lmm1 = _Q_l1;
	  auto _Q_lm = phase * (Tp{2} * x * _Q_lmm1 / root
				+ Tp(l + 1) * Tp(l) * _Q_l0);
	  for (unsigned int k = 3; k <= m; ++k)
	    {
	      _Q_lmm2 = _Q_lmm1;
	      _Q_lmm1 = _Q_lm;
	      _Q_lm = phase * (Tp(2 * (k - 1)) * x * _Q_lmm1 / root
			      - Tp(l + k - 1) * Tp(l - k + 2) * _Q_l0);
	    }

	  return {l, m, x, _Q_lm, _Q_lmm1, _Q_lmm2, phase};
	}
      else
	{ // FIXME!
	  return {l, m, x, Tp{0}, Tp{0}, Tp{0}, phase};
	}
    }

  /**
   * @brief  Return the spherical associated Legendre function.
   *
   * The spherical associated Legendre function of @f$ l @f$, @f$ m @f$,
   * and @f$ \theta @f$ is defined as @f$ Y_l^m(\theta,0) @f$ where
   * @f[
   * 	Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   * 				    \frac{(l-m)!}{(l+m)!}]
   * 		       P_l^m(\cos\theta) \exp^{im\phi}
   * @f]
   * is the spherical harmonic function and @f$ P_l^m(x) @f$ is the
   * associated Legendre function.
   *
   * This function differs from the associated Legendre function by
   * argument (@f$x = \cos(\theta)@f$) and by a normalization factor
   * but this factor is rather large for large @f$ l @f$ and @f$ m @f$
   * and so this function is stable for larger differences of @f$ l @f$
   * and @f$ m @f$.
   * @note Unlike the case for assoc_legendre_p the Condon-Shortley
   *       phase factor @f$ (-1)^m @f$ is present here.
   * @note @f$ Y_l^m(\theta) = 0 @f$ if @f$ m > l @f$.
   *
   * @param  l  The degree of the spherical associated Legendre function.
   * 		@f$ l >= 0 @f$.
   * @param  m  The order of the spherical associated Legendre function.
   * @param  theta  The radian polar angle argument
   * 		    of the spherical associated Legendre function.
   */
  template<typename Tp>
    Tp
    sph_legendre(unsigned int l, unsigned int m, Tp theta)
    {
      if (std::isnan(theta))
	return emsr::quiet_NaN(theta);

      const auto x = std::cos(theta);

      if (m > l)
	return Tp{0};
      else if (m == 0)
	{
	  auto _P_l = legendre_p(l, x).P_l;
	  Tp fact = std::sqrt(Tp(2 * l + 1)
		     / (Tp{4} * emsr::pi_v<Tp>));
	  _P_l *= fact;
	  return _P_l;
	}
      else if (x == Tp{1} || x == -Tp{1})
	return Tp{0}; // m > 0 here
      else
	{
	  // m > 0 and |x| < 1 here

	  // Starting value for recursion.
	  // Y_m^m(x) = \sqrt{ (2m+1)/(4pi m) \Gamma(m+1/2)/\Gamma(m) }
	  //           (-1)^m (1-x^2)^(m/2) / \pi^(1/4)
	  const auto sgn = (m % 2 == 1 ? -Tp{1} : Tp{1});
	  const auto _Y_mp1m_factor = x * std::sqrt(Tp(2 * m + 3));
	  const auto lncirc = std::log1p(-x * x);
	  // Gamma(m+1/2) / Gamma(m)
	  const auto lnpoch = log_gamma(Tp(m + 0.5L))
			      - log_gamma(Tp(m));
	  const auto lnpre_val =
		     -Tp{0.25L} * emsr::lnpi_v<Tp>
		     + Tp{0.5L} * (lnpoch + m * lncirc);
	  const auto sr = std::sqrt((Tp{2} + Tp{1} / m)
			  / (Tp{4} * emsr::pi_v<Tp>));
	  auto _Y_mm = sgn * sr * std::exp(lnpre_val);
	  auto _Y_mp1m = _Y_mp1m_factor * _Y_mm;

	  if (l == m)
	    return _Y_mm;
	  else if (l == m + 1)
	    return _Y_mp1m;
	  else
	    {
	      auto _Y_lm = Tp{0};

	      // Compute Y_l^m, l > m+1, upward recursion on l.
	      for (auto ll = m + 2; ll <= l; ++ll)
		{
		  const auto rat1 = Tp(ll - m) / Tp(ll + m);
		  const auto rat2 = Tp(ll - m - 1) / Tp(ll + m - 1);
		  const auto fact1 = std::sqrt(rat1 * Tp(2 * ll + 1)
						       * Tp(2 * ll - 1));
		  const auto fact2 = std::sqrt(rat1 * rat2
						 * Tp(2 * ll + 1)
						 / Tp(2 * ll - 3));
		  _Y_lm = (x * _Y_mp1m * fact1
			 - Tp(ll + m - 1) * _Y_mm * fact2)
			 / Tp(ll - m);
		  _Y_mm = _Y_mp1m;
		  _Y_mp1m = _Y_lm;
		}

	      return _Y_lm;
	    }
	}
    }


  /**
   * @brief  Return the spherical harmonic function.
   *
   * The spherical harmonic function of @f$ l @f$, @f$ m @f$,
   * and @f$ \theta @f$, @f$ \phi @f$ is defined by:
   * @f[
   * 	Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   * 				    \frac{(l-m)!}{(l+m)!}]
   * 		       P_l^{|m|}(\cos\theta) \exp^{im\phi}
   * @f]
   * @note @f$ Y_l^m(\theta,\phi) = 0 @f$ if @f$ |m| > l @f$.
   *
   * @param  l  The degree of the spherical harmonic function.
   * 		@f$ l >= 0 @f$.
   * @param  m  The order of the spherical harmonic function.
   * @param  theta  The radian polar angle argument
   * 		    of the spherical harmonic function.
   * @param  phi    The radian azimuthal angle argument
   * 		    of the spherical harmonic function.
   */
  template<typename Tp>
    std::complex<Tp>
    sph_harmonic(unsigned int l, int m, Tp theta, Tp phi)
    {
      const auto s_NaN = emsr::quiet_NaN(theta);
      if (std::isnan(theta) || std::isnan(phi))
	return std::complex<Tp>{s_NaN, s_NaN};
      else if (std::abs(m) > l)
	return std::complex<Tp>{0, 0};
      else
	return sph_legendre(l, std::abs(m), theta)
		* std::polar(Tp{1}, Tp(m) * phi);
    }


  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    legendre_zeros(unsigned int l, Tp proto = Tp{})
    {
      const auto s_eps = emsr::epsilon(proto);
      const auto s_pi = emsr::pi_v<Tp>;
      const unsigned int s_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<Tp>> pt(l);

      auto m = l / 2;

      // Treat the central zero for odd degree specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large degree.
      if (l & 1)
	{
	  const auto lm = l - 1;
	  const auto mm = lm / 2;
	  auto _Am = Tp{1};
	  for (auto m = 1u; m <= mm; ++m)
	    _Am *= -Tp(2 * m - 1) / Tp(2 * m);
	  auto Plm1 = _Am;
	  auto Ppl = l * Plm1;
	  pt[m].point = Tp{0};
	  pt[m].weight = Tp{2} / Ppl / Ppl;
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  // Clever approximation of root.
	  auto z = std::cos(s_pi * (i - Tp{1} / Tp{4})
				    / (l + Tp{1} / Tp{2}));
	  auto z1 = z;
	  auto w = Tp{0};
	  for (auto its = 0u; its < s_maxit; ++its)
	    {
	      // Compute P, P1, and P2 the Legendre polynomials of degree
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute Pp the derivative of the Legendre polynomial
	      // of degree l.
	      auto P1 = Tp{0};
	      auto P = Tp{1};
	      for  (auto k = 1u; k <= l; ++k)
		{
		  auto P2 = P1;
		  P1 = P;
		  // Recursion for Legendre polynomials.
		  P = ((Tp{2} * k - Tp{1}) * z * P1
		      - (k - Tp{1}) * P2) / k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto Pp = l * (z * P - P1) / (z * z - Tp{1});
	      z1 = z;
	      // Converge on root by Newton's method.
	      z = z1 - P / Pp;
	      if (std::abs(z - z1) < s_eps)
		{
		  w = Tp{2} / ((Tp{1} - z * z) * Pp * Pp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("legendre_zeros: Too many iterations");
	    }

	  pt[i - 1].point = -z;
	  pt[l - i].point = z;
	  pt[i - 1].weight = w;
	  pt[l - i].weight = w;
	}

      return pt;
    }

} // namespace detail
} // namespace emsr

#endif // SF_LEGENDRE_TCC
