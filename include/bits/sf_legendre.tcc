// Special functions -*- C++ -*-

// Copyright (C) 2006-2018 Free Software Foundation, Inc.
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
// (1) Handbook of Mathematical Functions,
//     ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 8, pp. 331-341
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 252-254

#ifndef _GLIBCXX_BITS_SF_LEGENDRE_TCC
#define _GLIBCXX_BITS_SF_LEGENDRE_TCC 1

#pragma GCC system_header

#include <complex>
#include <ext/math_const.h>
#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

// Implementation-space details.
namespace __detail
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
   * @param  __l  The degree of the Legendre polynomial.  @f$ l >= 0 @f$.
   * @param  __x  The argument of the Legendre polynomial.
   */
  template<typename _Tp>
    __gnu_cxx::__legendre_p_t<_Tp>
    __legendre_p(unsigned int __l, _Tp __x)
    {
      using __ret_t = __gnu_cxx::__legendre_p_t<_Tp>;

      const auto __lge1 = __l >= 1 ? _Tp{+1} : _Tp{0};
      const auto __lge2 = __l >= 2 ? _Tp{+1} : _Tp{0};
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);

      if (std::isnan(__x))
	return {__l, _S_NaN, _S_NaN, _S_NaN, _S_NaN};
      else if (__x == _Tp{+1})
	return {__l, __x, _Tp{+1}, __lge1, __lge2};
      else if (__x == _Tp{-1})
	return __l % 2 == 1
		? __ret_t{__l, __x, _Tp{-1}, +__lge1, -__lge2}
		: __ret_t{__l, __x, _Tp{+1}, -__lge1, +__lge2};
      else
	{
	  auto _P_lm2 = _Tp{1};
	  if (__l == 0)
	    return {__l, __x, _P_lm2, _Tp{0}, _Tp{0}};

	  auto _P_lm1 = __x;
	  if (__l == 1)
	    return {__l, __x, _P_lm1, _P_lm2, _Tp{0}};

	  auto _P_l = _Tp{2} * __x * _P_lm1 - _P_lm2
		    - (__x * _P_lm1 - _P_lm2) / _Tp{2};
	  for (unsigned int __ll = 3; __ll <= __l; ++__ll)
	    {
	      _P_lm2 = _P_lm1;
	      _P_lm1 = _P_l;
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      _P_l = _Tp{2} * __x * _P_lm1 - _P_lm2
		    - (__x * _P_lm1 - _P_lm2) / _Tp(__ll);
	    }
	  // Recursion for the derivative of The Legendre polynomial.
	  //auto __Pp_l = __l * (__z * _P_l - _P_lm1) / (__z * __z - _Tp{1});

	  return {__l, __x, _P_l, _P_lm1, _P_lm2};
	}
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
   * @param __l The degree of the Legendre function.  @f$l >= 0@f$.
   * @param __x The argument of the Legendre function.  @f$|x| <= 1@f$.
   */
  template<typename _Tp>
    _Tp
    __legendre_q(unsigned int __l, _Tp __x)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      if ((__x < -_Tp{1}) || (__x > +_Tp{1}))
	std::__throw_domain_error(__N("__legendre_q: argument out of range"));
      else if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (std::abs(__x - _Tp{1}) < _S_eps)
	return _S_inf;
      else if (std::abs(__x + _Tp{1}) < _S_eps)
	return (__l & 1 ? +1 : -1) * _S_inf;
      else
	{
	  auto _Q_lm2 = _Tp{0.5L} * std::log((_Tp{1} + __x) / (_Tp{1} - __x));
	  if (__l == 0)
	    return _Q_lm2;
	  auto _Q_lm1 = __x * _Q_lm2 - _Tp{1};
	  if (__l == 1)
	    return _Q_lm1;
	  auto _Q_l = _Tp{2} * __x * _Q_lm1 - _Q_lm2
		    - (__x * _Q_lm1 - _Q_lm2) / _Tp{2};
	  for (unsigned int __ll = 3; __ll <= __l; ++__ll)
	    {
	      _Q_lm2 = _Q_lm1;
	      _Q_lm1 = _Q_l;
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      _Q_l = _Tp{2} * __x * _Q_lm1 - _Q_lm2
		    - (__x * _Q_lm1 - _Q_lm2) / _Tp(__ll);
	    }

	  return _Q_l;
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
   * by default.
   *
   * @param  __l  The degree of the associated Legendre function.
   * 		@f$ l >= 0 @f$.
   * @param  __m  The order of the associated Legendre function.
   * 		@f$ m <= l @f$.
   * @param  __x  The argument of the associated Legendre function.
   * @param  __phase  The phase of the associated Legendre function.
   *                  Use -1 for the Condon-Shortley phase convention.
   */
  template<typename _Tp>
    _Tp
    __assoc_legendre_p(unsigned int __l, unsigned int __m, _Tp __x,
		       _Tp __phase = _Tp{+1})
    {
      if (__m > __l)
	std::__throw_domain_error(__N("__assoc_legendre_p: "
				      "degree out of range"));
      else if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__m == 0)
	return __legendre_p(__l, __x).__P_l;
      else
	{
	  _Tp _P_mm = _Tp{1};
	  if (__m > 0)
	    {
	      // Two square roots seem more accurate more of the time
	      // than just one.
	      auto __root = std::sqrt(_Tp{1} - __x) * std::sqrt(_Tp{1} + __x);
	      auto __fact = _Tp{1};
	      for (unsigned int __i = 1; __i <= __m; ++__i)
		{
		  // N. B. Condon-Shortley would use __phase = -1 here.
		  _P_mm *= __phase * __fact * __root;
		  __fact += _Tp{2};
		}
	    }
	  if (__l == __m)
	    return _P_mm;

	  _Tp _P_mp1m = _Tp(2 * __m + 1) * __x * _P_mm;
	  if (__l == __m + 1)
	    return _P_mp1m;

	  auto _P_lm2m = _P_mm;
	  auto _P_lm1m = _P_mp1m;
	  auto _P_lm = (_Tp(2 * __m + 3) * __x * _P_lm1m
		     - _Tp(2 * __m + 1) * _P_lm2m) / _Tp{2};
	  for (unsigned int __j = __m + 3; __j <= __l; ++__j)
	    {
	      _P_lm2m = _P_lm1m;
	      _P_lm1m = _P_lm;
	      _P_lm = (_Tp(2 * __j - 1) * __x * _P_lm1m
		      - _Tp(__j + __m - 1) * _P_lm2m) / _Tp(__j - __m);
	    }

	  return _P_lm;
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
   * @note Unlike the case for __assoc_legendre_p the Condon-Shortley
   * phase factor @f$ (-1)^m @f$ is present here.
   *
   * @param  __l  The degree of the spherical associated Legendre function.
   * 		@f$ l >= 0 @f$.
   * @param  __m  The order of the spherical associated Legendre function.
   * 		@f$ m <= l @f$.
   * @param  __theta  The radian polar angle argument
   * 		    of the spherical associated Legendre function.
   */
  template<typename _Tp>
    _Tp
    __sph_legendre(unsigned int __l, unsigned int __m, _Tp __theta)
    {
      if (std::isnan(__theta))
	return __gnu_cxx::__quiet_NaN(__theta);

      const auto __x = std::cos(__theta);

      if (__l < __m)
	std::__throw_domain_error(__N("__sph_legendre: bad argument"));
      else if (__m == 0)
	{
	  auto _P_l = __legendre_p(__l, __x).__P_l;
	  _Tp __fact = std::sqrt(_Tp(2 * __l + 1)
		     / (_Tp{4} * __gnu_cxx::__const_pi(__theta)));
	  _P_l *= __fact;
	  return _P_l;
	}
      else if (__x == _Tp{1} || __x == -_Tp{1})
	return _Tp{0}; // m > 0 here
      else
	{
	  // m > 0 and |x| < 1 here

	  // Starting value for recursion.
	  // Y_m^m(x) = \sqrt{ (2m+1)/(4pi m) \Gamma(m+1/2)/\Gamma(m) }
	  //           (-1)^m (1-x^2)^(m/2) / \pi^(1/4)
	  const auto __sgn = (__m % 2 == 1 ? -_Tp{1} : _Tp{1});
	  const auto _Y_mp1m_factor = __x * std::sqrt(_Tp(2 * __m + 3));
	  const auto __lncirc = std::log1p(-__x * __x);
	  // Gamma(m+1/2) / Gamma(m)
	  const auto __lnpoch = __log_gamma(_Tp(__m + 0.5L))
			      - __log_gamma(_Tp(__m));
	  const auto __lnpre_val =
		     -_Tp{0.25L} * __gnu_cxx::__const_ln_pi(__theta)
		     + _Tp{0.5L} * (__lnpoch + __m * __lncirc);
	  const auto __sr = std::sqrt((_Tp{2} + _Tp{1} / __m)
			  / (_Tp{4} * __gnu_cxx::__const_pi(__theta)));
	  auto _Y_mm = __sgn * __sr * std::exp(__lnpre_val);
	  auto _Y_mp1m = _Y_mp1m_factor * _Y_mm;

	  if (__l == __m)
	    return _Y_mm;
	  else if (__l == __m + 1)
	    return _Y_mp1m;
	  else
	    {
	      auto _Y_lm = _Tp{0};

	      // Compute Y_l^m, l > m+1, upward recursion on l.
	      for ( int __ll = __m + 2; __ll <= __l; ++__ll)
		{
		  const auto __rat1 = _Tp(__ll - __m) / _Tp(__ll + __m);
		  const auto __rat2 = _Tp(__ll - __m - 1) / _Tp(__ll + __m - 1);
		  const auto __fact1 = std::sqrt(__rat1 * _Tp(2 * __ll + 1)
						       * _Tp(2 * __ll - 1));
		  const auto __fact2 = std::sqrt(__rat1 * __rat2
						 * _Tp(2 * __ll + 1)
						 / _Tp(2 * __ll - 3));
		  _Y_lm = (__x * _Y_mp1m * __fact1
			 - _Tp(__ll + __m - 1) * _Y_mm * __fact2)
			 / _Tp(__ll - __m);
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
   *
   * @param  __l  The degree of the spherical harmonic function.
   * 		@f$ l >= 0 @f$.
   * @param  __m  The order of the spherical harmonic function.
   * 		@f$ m <= l @f$.
   * @param  __theta  The radian polar angle argument
   * 		    of the spherical harmonic function.
   * @param  __phi    The radian azimuthal angle argument
   * 		    of the spherical harmonic function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sph_harmonic(unsigned int __l, int __m, _Tp __theta, _Tp __phi)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__theta);
      if (std::isnan(__theta) || std::isnan(__phi))
	return std::complex<_Tp>{_S_NaN, _S_NaN};

      return __sph_legendre(__l, std::abs(__m), __theta)
	   * std::polar(_Tp{1}, _Tp(__m) * __phi);
    }


  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __legendre_zeros(unsigned int __l, _Tp proto = _Tp{})
    {
      const auto _S_eps = __gnu_cxx::__epsilon(proto);
      const auto _S_pi = __gnu_cxx::__const_pi(proto);
      const unsigned int _S_maxit = 1000u;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__l);

      auto __m = __l / 2;

      // Treat the central zero for odd degree specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large degree.
      if (__l & 1)
	{
	  const auto __lm = __l - 1;
	  const auto __mm = __lm / 2;
	  auto _Am = _Tp{1};
	  for (auto __m = 1u; __m <= __mm; ++__m)
	    _Am *= -_Tp(2 * __m - 1) / _Tp(2 * __m);
	  auto __Plm1 = _Am;
	  auto __Ppl = __l * __Plm1;
	  __pt[__m].__point = _Tp{0};
	  __pt[__m].__weight = _Tp{2} / __Ppl / __Ppl;
	}

      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  // Clever approximation of root.
	  auto __z = std::cos(_S_pi * (__i - _Tp{1} / _Tp{4})
				    / (__l + _Tp{1} / _Tp{2}));
	  auto __z1 = __z;
	  auto __w = _Tp{0};
	  for (auto __its = 0u; __its < _S_maxit; ++__its)
	    {
	      // Compute __P, __P1, and __P2 the Legendre polynomials of degree
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute __Pp the derivative of the Legendre polynomial
	      // of degree l.
	      auto __P1 = _Tp{0};
	      auto __P = _Tp{1};
	      for  (auto __k = 1u; __k <= __l; ++__k)
		{
		  auto __P2 = __P1;
		  __P1 = __P;
		  // Recursion for Legendre polynomials.
		  __P = ((_Tp{2} * __k - _Tp{1}) * __z * __P1
		      - (__k - _Tp{1}) * __P2) / __k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto __Pp = __l * (__z * __P - __P1) / (__z * __z - _Tp{1});
	      __z1 = __z;
	      // Converge on root by Newton's method.
	      __z = __z1 - __P / __Pp;
	      if (std::abs(__z - __z1) < _S_eps)
		{
		  __w = _Tp{2} / ((_Tp{1} - __z * __z) * __Pp * __Pp);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__legendre_zeros: "
					 "Too many iterations");
	    }

	  __pt[__i - 1].__point = -__z;
	  __pt[__l - __i].__point = __z;
	  __pt[__i - 1].__weight = __w;
	  __pt[__l - __i].__weight = __w;
	}

      return __pt;
    }
} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_LEGENDRE_TCC