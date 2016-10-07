// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
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

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *   @brief  Return the Legendre polynomial by upward recursion
   *           on order @f$ l @f$.
   *
   *   The Legendre function of order @f$ l @f$ and argument @f$ x @f$,
   *   @f$ P_l(x) @f$, is defined by:
   *   @f[
   *     P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   *   @f]
   *
   *   @param  __l  The order of the Legendre polynomial.  @f$l >= 0@f$.
   *   @param  __x  The argument of the Legendre polynomial.  @f$|x| <= 1@f$.
   */
  template<typename _Tp>
    _Tp
    __poly_legendre_p(unsigned int __l, _Tp __x)
    {
      if ((__x < -_Tp{1}) || (__x > +_Tp{1}))
	std::__throw_domain_error(__N("__poly_legendre_p: argument out of range"));
      else if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__x == +_Tp{1})
	return +_Tp{1};
      else if (__x == -_Tp{1})
	return (__l % 2 == 1 ? -_Tp{1} : +_Tp{1});
      else
	{
	  auto _P_lm2 = _Tp{1};
	  if (__l == 0)
	    return _P_lm2;

	  auto _P_lm1 = __x;
	  if (__l == 1)
	    return _P_lm1;

	  auto _P_l = _Tp{0};
	  for (unsigned int __ll = 2; __ll <= __l; ++__ll)
	    {
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      _P_l = _Tp{2} * __x * _P_lm1 - _P_lm2
		    - (__x * _P_lm1 - _P_lm2) / _Tp(__ll);
	      _P_lm2 = _P_lm1;
	      _P_lm1 = _P_l;
	    }

	  return _P_l;
	}
    }

  /**
   * @brief Return the Legendre function of the second kind
   *        by upward recursion on order @f$ l @f$.
   *
   * The Legendre function of the second kind of order @f$ l @f$
   * and argument @f$ x @f$, @f$ Q_l(x) @f$, is defined by:
   * @f[
   *   Q_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   *
   * @param __l The order of the Legendre function.  @f$l >= 0@f$.
   * @param __x The argument of the Legendre function.  @f$|x| <= 1@f$.
   */
  template<typename _Tp>
    _Tp
    __legendre_q(unsigned int __l, _Tp __x)
    {
      if ((__x < -_Tp{1}) || (__x > +_Tp{1}))
	std::__throw_domain_error(__N("__legendre_q: argument out of range"));
      else if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__x == +_Tp{1})
	return +_Tp{1};
      else if (__x == -_Tp{1})
	return (__l % 2 == 1 ? -_Tp{1} : +_Tp{1});
      else
	{
	  auto _Q_lm2 = _Tp{0.5L} * std::log((_Tp{1} + __x) / (_Tp{1} - __x));
	  if (__l == 0)
	    return _Q_lm2;
	  auto _Q_lm1 = __x * _Q_lm2 - _Tp{1};
	  if (__l == 1)
	    return _Q_lm1;
	  auto _Q_l = _Tp{0};
	  for (unsigned int __ll = 2; __ll <= __l; ++__ll)
	    {
	      // This arrangement is supposed to be better for roundoff
	      // protection, Arfken, 2nd Ed, Eq 12.17a.
	      _Q_l = _Tp{2} * __x * _Q_lm1 - _Q_lm2
		    - (__x * _Q_lm1 - _Q_lm2) / _Tp(__ll);
	      _Q_lm2 = _Q_lm1;
	      _Q_lm1 = _Q_l;
	    }

	  return _Q_l;
	}
    }

  /**
   *   @brief  Return the associated Legendre function by recursion
   *           on @f$ l @f$ and downward recursion on m.
   *
   *   The associated Legendre function is derived from the Legendre function
   *   @f$ P_l(x) @f$ by the Rodrigues formula:
   *   @f[
   *     P_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
   *   @f]
   *
   *   @param  __l  The order of the associated Legendre function.
   *              @f$ l >= 0 @f$.
   *   @param  __m  The order of the associated Legendre function.
   *              @f$ m <= l @f$.
   *   @param  __x  The argument of the associated Legendre function.
   *              @f$ |x| <= 1 @f$.
   */
  template<typename _Tp>
    _Tp
    __assoc_legendre_p(unsigned int __l, unsigned int __m, _Tp __x)
    {
      if (__x < -_Tp{1} || __x > +_Tp{1})
	std::__throw_domain_error(__N("__assoc_legendre_p: "
				      "argument out of range"));
      else if (__m > __l)
	std::__throw_domain_error(__N("__assoc_legendre_p: "
				      "degree out of range"));
      else if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__m == 0)
	return __poly_legendre_p(__l, __x);
      else
	{
	  _Tp _P_mm = _Tp{1};
	  if (__m > 0)
	    {
	      // Two square roots seem more accurate more of the time
	      // than just one.
	      _Tp __root = std::sqrt(_Tp{1} - __x) * std::sqrt(_Tp{1} + __x);
	      _Tp __fact = _Tp{1};
	      for (unsigned int __i = 1; __i <= __m; ++__i)
		{
		  _P_mm *= -__fact * __root;
		  __fact += _Tp{2};
		}
	    }
	  if (__l == __m)
	    return _P_mm;

	  _Tp _P_mp1m = _Tp(2 * __m + 1) * __x * _P_mm;
	  if (__l == __m + 1)
	    return _P_mp1m;

	  _Tp _P_lm2m = _P_mm;
	  _Tp _P_lm1m = _P_mp1m;
	  _Tp _P_lm = _Tp{0};
	  for (unsigned int __j = __m + 2; __j <= __l; ++__j)
	    {
	      _P_lm = (_Tp(2 * __j - 1) * __x * _P_lm1m
		      - _Tp(__j + __m - 1) * _P_lm2m) / _Tp(__j - __m);
	      _P_lm2m = _P_lm1m;
	      _P_lm1m = _P_lm;
	    }

	  return _P_lm;
	}
    }


  /**
   *   @brief  Return the spherical associated Legendre function.
   *
   *   The spherical associated Legendre function of @f$ l @f$, @f$ m @f$,
   *   and @f$ \theta @f$ is defined as @f$ Y_l^m(\theta,0) @f$ where
   *   @f[
   *      Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   *                                  \frac{(l-m)!}{(l+m)!}]
   *                     P_l^m(\cos\theta) \exp^{im\phi}
   *   @f]
   *   is the spherical harmonic function and @f$ P_l^m(x) @f$ is the
   *   associated Legendre function.
   *
   *   This function differs from the associated Legendre function by
   *   argument (@f$x = \cos(\theta)@f$) and by a normalization factor
   *   but this factor is rather large for large @f$ l @f$ and @f$ m @f$
   *   and so this function is stable for larger differences of @f$ l @f$
   *   and @f$ m @f$.
   *
   *   @param  __l  The order of the spherical associated Legendre function.
   *              @f$ l >= 0 @f$.
   *   @param  __m  The order of the spherical associated Legendre function.
   *              @f$ m <= l @f$.
   *   @param  __theta  The radian polar angle argument
   *                  of the spherical associated Legendre function.
   */
  template<typename _Tp>
    _Tp
    __sph_legendre(unsigned int __l, unsigned int __m, _Tp __theta)
    {
      if (__isnan(__theta))
	return __gnu_cxx::__quiet_NaN<_Tp>();

      const auto __x = std::cos(__theta);

      if (__l < __m)
	std::__throw_domain_error(__N("__sph_legendre: bad argument"));
      else if (__m == 0)
	{
	  _Tp _P_l = __poly_legendre_p(__l, __x);
	  _Tp __fact = std::sqrt(_Tp(2 * __l + 1)
		     / (_Tp{4} * __gnu_cxx::__math_constants<_Tp>::__pi));
	  _P_l *= __fact;
	  return _P_l;
	}
      else if (__x == _Tp{1} || __x == -_Tp{1})
	return _Tp{0}; // m > 0 here
      else
	{
	  // m > 0 and |x| < 1 here

	  // Starting value for recursion.
	  // Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) )
	  //           (-1)^m (1-x^2)^(m/2) / pi^(1/4)
	  const auto __sgn = (__m % 2 == 1 ? -_Tp{1} : _Tp{1});
	  const auto _Y_mp1m_factor = __x * std::sqrt(_Tp(2 * __m + 3));
	  const auto __lncirc = std::log1p(-__x * __x);
	  // Gamma(m+1/2) / Gamma(m)
	  const auto __lnpoch = __log_gamma(_Tp(__m + 0.5L))
			      - __log_gamma(_Tp(__m));
	  const auto __lnpre_val =
		     -_Tp{0.25L} * __gnu_cxx::__math_constants<_Tp>::__ln_pi
		     + _Tp{0.5L} * (__lnpoch + __m * __lncirc);
	  _Tp __sr = std::sqrt((_Tp{2} + _Tp{1} / __m)
		   / (_Tp{4} * __gnu_cxx::__math_constants<_Tp>::__pi));
	  _Tp _Y_mm = __sgn * __sr * std::exp(__lnpre_val);
	  _Tp _Y_mp1m = _Y_mp1m_factor * _Y_mm;

	  if (__l == __m)
	    {
	      return _Y_mm;
	    }
	  else if (__l == __m + 1)
	    {
	      return _Y_mp1m;
	    }
	  else
	    {
	      _Tp _Y_lm = _Tp{0};

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
   *   @brief  Return the spherical harmonic function.
   *
   *   The spherical harmonic function of @f$ l @f$, @f$ m @f$,
   *   and @f$ \theta @f$, @f$ \phi @f$ is defined by:
   *   @f[
   *      Y_l^m(\theta,\phi) = (-1)^m[\frac{(2l+1)}{4\pi}
   *                                  \frac{(l-m)!}{(l+m)!}]
   *                     P_l^{|m|}(\cos\theta) \exp^{im\phi}
   *   @f]
   *
   *   @param  __l  The order of the spherical harmonic function.
   *              @f$ l >= 0 @f$.
   *   @param  __m  The order of the spherical harmonic function.
   *              @f$ m <= l @f$.
   *   @param  __theta  The radian polar angle argument
   *                  of the spherical harmonic function.
   *   @param  __phi    The radian azimuthal angle argument
   *                  of the spherical harmonic function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sph_harmonic(unsigned int __l, int __m, _Tp __theta, _Tp __phi)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();
      if (__isnan(__theta) || __isnan(__phi))
	return std::complex<_Tp>{_S_NaN, _S_NaN};

      return __sph_legendre(__l, std::abs(__m), __theta)
	   * std::polar(_Tp{1}, _Tp(__m) * __phi);
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_LEGENDRE_TCC
