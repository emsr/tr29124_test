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

#ifndef _GLIBCXX_BITS_SF_BESSEL_TCC
#define _GLIBCXX_BITS_SF_BESSEL_TCC 1

#pragma GCC system_header

#include <complex>
#include <utility> // For exchange
#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief This routine computes the asymptotic cylindrical Bessel
   * 	    and Neumann functions of order nu: @f$ J_{\nu}(x) @f$,
   * 	    @f$ N_{\nu}(x) @f$.  Use this for @f$ x >> nu^2 + 1 @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @param[out]  _Jnu  The Bessel function of the first kind.
   * @param[out]  _Nnu  The Neumann function (Bessel function of the second kind).
   * @param[out]  _Jpnu  The Bessel function of the first kind.
   * @param[out]  _Npnu  The Neumann function (Bessel function of the second kind).
   */
  template<typename _Tp>
    void
    __cyl_bessel_jn_asymp(_Tp __nu, _Tp __x,
			  _Tp& _Jnu, _Tp& _Nnu,
			  _Tp& _Jpnu, _Tp& _Npnu)
    {
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<long double>::__pi_half;
      const auto __2nu = _Tp{2} * __nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __8x = _Tp{8} * __x;
      auto __k = 0;
      auto __bk_xk = _Tp{1};
      auto _Rsum = __bk_xk;
      auto __ak_xk = _Tp{1};
      auto _Psum = __ak_xk;
      ++__k;
      auto __2km1 = 1;
      __bk_xk *= (__4nu2 + __2km1 * (__2km1 + 2)) / __8x;
      auto _Ssum = __bk_xk;
      __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / __8x;
      auto _Qsum = __ak_xk;
      do
	{
	  ++__k;
	  __2km1 += 2;
	  __bk_xk = -(__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Rsum += __bk_xk;
	  __ak_xk *= -(__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
	  _Psum += __ak_xk;
	  auto __convP = std::abs(__ak_xk) < _S_eps * std::abs(_Psum);

	  ++__k;
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Ssum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
	  _Qsum += __ak_xk;
	  auto __convQ = std::abs(__ak_xk) < _S_eps * std::abs(_Qsum);

	  if (__convP && __convQ && __k > (__nu / _Tp{2}))
	    break;
	}
      while (__k < _Tp{100} * __nu);

      auto __omega = __x - (__nu + 0.5L) * _S_pi_2;
      auto __c = std::cos(__omega);
      auto __s = std::sin(__omega);

      auto __coef = std::sqrt(_Tp{2} / (_S_pi * __x));
      _Jnu = __coef * (__c * _Psum - __s * _Qsum);
      _Nnu = __coef * (__s * _Psum + __c * _Qsum);
      _Jpnu = -__coef * (__s * _Rsum + __c * _Ssum);
      _Npnu =  __coef * (__c * _Rsum - __s * _Ssum);

      return;
    }

  /**
   * @brief Compute the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 		  [\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   * 		  [\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are exquisite.
   *
   * @param __mu     The input parameter of the gamma functions.
   * @param[out] __gam1   The output function @f$ \Gamma_1(\mu) @f$
   * @param[out] __gam2   The output function @f$ \Gamma_2(\mu) @f$
   * @param[out] __gampl  The output function @f$ \Gamma(1 + \mu) @f$
   * @param[out] __gammi  The output function @f$ \Gamma(1 - \mu) @f$
   */
  template<typename _Tp>
    void
    __gamma_temme(_Tp __mu,
		  _Tp& __gam1, _Tp& __gam2, _Tp& __gampl, _Tp& __gammi)
    {
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_gamma_E = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      __gampl = _Tp{1} / __gamma(_Tp{1} + __mu);
      __gammi = _Tp{1} / __gamma(_Tp{1} - __mu);

      if (std::abs(__mu) < _S_eps)
	__gam1 = -_S_gamma_E;
      else
	__gam1 = (__gammi - __gampl) / (_Tp{2} * __mu);

      __gam2 = (__gammi + __gampl) / _Tp{2};

      return;
    }

  /**
   * @brief  Compute the Bessel @f$ J_\nu(x) @f$ and Neumann
   * 	     @f$ N_\nu(x) @f$ functions and their first derivatives
   * 	     @f$ J'_\nu(x) @f$ and @f$ N'_\nu(x) @f$ respectively.
   * 	     These four functions are computed together for numerical
   * 	     stability.
   *
   * @param __nu The order of the Bessel functions.
   * @param __x  The argument of the Bessel functions.
   * @param[out] _Jnu The output Bessel function of the first kind.
   * @param[out] _Nnu The output Neumann function (Bessel function of the second kind).
   * @param[out] _Jpnu The output derivative of the Bessel function of the first kind.
   * @param[out] _Npnu The output derivative of the Neumann function.
   */
  template<typename _Tp>
    void
    __cyl_bessel_jn_steed(_Tp __nu, _Tp __x,
			  _Tp& _Jnu, _Tp& _Nnu, _Tp& _Jpnu, _Tp& _Npnu)
    {
      constexpr auto _S_inf = __gnu_cxx::__infinity<_Tp>();
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      // When the multiplier is N i.e.
      // fp_min = N * min()
      // Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
      //const _Tp _S_fp_min = _Tp{20} * __gnu_cxx::__min<_Tp>();
      constexpr int _S_max_iter = 15000;
      constexpr auto _S_x_min = _Tp{2};
      const auto _S_fp_min = __gnu_cxx::__sqrt_min<_Tp>();

      const int __n = (__x < _S_x_min
		    ? std::nearbyint(__nu)
		    : std::max(0,
			       static_cast<int>(__nu - __x + _Tp{1.5L})));

      const auto __mu = __nu - __n;
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      const auto _Wronski = __xi2 / _S_pi;
      int __isign = 1;
      _Tp __h = __nu * __xi;
      if (__h < _S_fp_min)
	__h = _S_fp_min;
      auto __b = __xi2 * __nu;
      auto __d = _Tp{0};
      auto __c = __h;
      int __i;
      for (__i = 1; __i <= _S_max_iter; ++__i)
	{
	  __b += __xi2;
	  __d = __b - __d;
	  if (std::abs(__d) < _S_fp_min)
	    __d = _S_fp_min;
	  __c = __b - _Tp{1} / __c;
	  if (std::abs(__c) < _S_fp_min)
	    __c = _S_fp_min;
	  __d = _Tp{1} / __d;
	  const auto __del = __c * __d;
	  __h *= __del;
	  if (__d < _Tp{0})
	    __isign = -__isign;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	}
      if (__i > _S_max_iter)
	{
	  // Don't throw with message "try asymptotic expansion" - Just do it!
	  __cyl_bessel_jn_asymp(__nu, __x, _Jnu, _Nnu, _Jpnu, _Npnu);
	  return;
	}
      auto _Jnul = __isign * _S_fp_min;
      auto _Jpnul = __h * _Jnul;
      auto _Jnul1 = _Jnul;
      auto _Jpnu1 = _Jpnul;
      auto __fact = __nu * __xi;
      for (int __l = __n; __l >= 1; --__l)
	{
	  const auto _Jnutemp = __fact * _Jnul + _Jpnul;
	  __fact -= __xi;
	  _Jpnul = __fact * _Jnutemp - _Jnul;
	  _Jnul = _Jnutemp;
	}
      if (_Jnul == _Tp{0})
	_Jnul = _S_eps;

      auto __f = _Jpnul / _Jnul;
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (__x < _S_x_min)
	{
	  const auto __x2 = __x / _Tp{2};
	  const auto __pimu = _S_pi * __mu;
	  const auto __fact = (std::abs(__pimu) < _S_eps
			    ? _Tp{1}
			    : __pimu / std::sin(__pimu));
	  auto __d = -std::log(__x2);
	  auto __e = __mu * __d;
	  const auto __fact2 = (std::abs(__e) < _S_eps
			     ? _Tp{1}
			     : std::sinh(__e) / __e);
	  _Tp __gam1, __gam2, __gampl, __gammi;
	  __gamma_temme(__mu, __gam1, __gam2, __gampl, __gammi);
	  auto __ff = (_Tp{2} / _S_pi) * __fact
		   * (__gam1 * std::cosh(__e) + __gam2 * __fact2 * __d);
	  __e = std::exp(__e);
	  auto __p = __e / (_S_pi * __gampl);
	  auto __q = _Tp{1} / (__e * _S_pi * __gammi);
	  const auto __pimu2 = __pimu / _Tp{2};
	  auto __fact3 = (std::abs(__pimu2) < _S_eps
		       ? _Tp{1} : std::sin(__pimu2) / __pimu2 );
	  auto __r = _S_pi * __pimu2 * __fact3 * __fact3;
	  auto __c = _Tp{1};
	  __d = -__x2 * __x2;
	  auto __sum = __ff + __r * __q;
	  auto __sum1 = __p;
	  int __i;
	  for (__i = 1; __i <= _S_max_iter; ++__i)
	    {
	      __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	      __c *= __d / _Tp(__i);
	      __p /= _Tp(__i) - __mu;
	      __q /= _Tp(__i) + __mu;
	      const auto __del = __c * (__ff + __r * __q);
	      __sum += __del;
	      const auto __del1 = __c * __p - _Tp(__i) * __del;
	      __sum1 += __del1;
	      if (std::abs(__del) < _S_eps * (_Tp{1} + std::abs(__sum)))
		break;
	    }
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__cyl_bessel_jn_steed: "
					   "Y-series failed to converge"));
	  _Nmu = -__sum;
	  _Nnu1 = -__sum1 * __xi2;
	  _Npmu = __mu * __xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - __f * _Nmu);
	}
      else
	{
	  auto __a = _Tp{0.25L} - __mu2;
	  auto __q = _Tp{1};
	  auto __p = -__xi / _Tp{2};
	  auto __br = _Tp{2} * __x;
	  auto __bi = _Tp{2};
	  auto __fact = __a * __xi / (__p * __p + __q * __q);
	  auto __cr = __br + __q * __fact;
	  auto __ci = __bi + __p * __fact;
	  auto __den = __br * __br + __bi * __bi;
	  auto __dr = __br / __den;
	  auto __di = -__bi / __den;
	  auto __dlr = __cr * __dr - __ci * __di;
	  auto __dli = __cr * __di + __ci * __dr;
	  auto __temp = __p * __dlr - __q * __dli;
	  __q = __p * __dli + __q * __dlr;
	  __p = __temp;
	  int __i;
	  for (__i = 2; __i <= _S_max_iter; ++__i)
	    {
	      __a += _Tp{2 * (__i - 1)};
	      __bi += _Tp{2};
	      __dr = __a * __dr + __br;
	      __di = __a * __di + __bi;
	      if (std::abs(__dr) + std::abs(__di) < _S_fp_min)
		__dr = _S_fp_min;
	      __fact = __a / (__cr * __cr + __ci * __ci);
	      __cr = __br + __cr * __fact;
	      __ci = __bi - __ci * __fact;
	      if (std::abs(__cr) + std::abs(__ci) < _S_fp_min)
		__cr = _S_fp_min;
	      __den = __dr * __dr + __di * __di;
	      __dr /= __den;
	      __di /= -__den;
	      __dlr = __cr * __dr - __ci * __di;
	      __dli = __cr * __di + __ci * __dr;
	      __temp = __p * __dlr - __q * __dli;
	      __q = __p * __dli + __q * __dlr;
	      __p = __temp;
	      if (std::abs(__dlr - _Tp{1}) + std::abs(__dli) < _S_eps)
		break;
	    }
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__cyl_bessel_jn_steed: "
					   "Lentz's method failed"));
	  const auto __gam = (__p - __f) / __q;
	  _Jmu = std::sqrt(_Wronski / ((__p - __f) * __gam + __q));
	  _Jmu = std::copysign(_Jmu, _Jnul);
	  _Nmu = __gam * _Jmu;
	  _Npmu = (__p + __q / __gam) * _Nmu;
	  _Nnu1 = __mu * __xi * _Nmu - _Npmu;
        }
      __fact = _Jmu / _Jnul;
      _Jnu = __fact * _Jnul1;
      _Jpnu = __fact * _Jpnu1;
      for (int __i = 1; __i <= __n; ++__i)
	_Nmu = std::exchange(_Nnu1, (__mu + __i) * __xi2 * _Nnu1 - _Nmu);
      _Nnu = _Nmu;
      _Npnu = __nu * __xi * _Nmu - _Nnu1;

      return;
    }


  /**
   * @brief This routine returns the cylindrical Bessel functions
   * 	    of order @f$ \nu @f$: @f$ J_{\nu} @f$ or @f$ I_{\nu} @f$
   * 	    by series expansion.
   *
   * The modified cylindrical Bessel function is:
   * @f[
   *  Z_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{\sigma^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   * where @f$ \sigma = +1 @f$ or@f$  -1 @f$ for
   * @f$ Z = I @f$ or @f$ J @f$ respectively.
   *
   * See Abramowitz & Stegun, 9.1.10
   * 	 Abramowitz & Stegun, 9.6.7
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Equation 9.1.10 p. 360 and Equation 9.6.10 p. 375
   *
   * @param  __nu  The order of the Bessel function.
   * @param  __x   The argument of the Bessel function.
   * @param  __sgn  The sign of the alternate terms
   * 		    -1 for the Bessel function of the first kind.
   * 		    +1 for the modified Bessel function of the first kind.
   * @param  __max_iter  The maximum number of iterations for sum.
   * @return  The output Bessel function.
   */
  template<typename _Tp>
    _Tp
    __cyl_bessel_ij_series(_Tp __nu, _Tp __x, _Tp __sgn,
			   unsigned int __max_iter)
    {
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      if (__x < _S_eps)
	{
	  if (__nu == _Tp{0})
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	{
	  const auto __x2 = __x / _Tp{2};

	  _Tp __fact = __nu * std::log(__x2);
	  __fact -= __log_gamma(__nu + _Tp{1});
	  __fact = std::exp(__fact);
	  const auto __xx4 = __sgn * __x2 * __x2;
	  _Tp _Jn = _Tp{1};
	  _Tp __term = _Tp{1};
	  for (unsigned int __i = 1; __i < __max_iter; ++__i)
	    {
	      __term *= __xx4 / (_Tp(__i) * (__nu + _Tp(__i)));
	      _Jn += __term;
	      if (std::abs(__term / _Jn) < _S_eps)
		break;
	    }

	  return __fact * _Jn;
	}
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   * of order @f$ \nu @f$ by various means.
   */
  template<typename _Tp>
    void
    __cyl_bessel_jn(_Tp __nu, _Tp __x,
		_Tp& _Jnu, _Tp& _Nnu, _Tp& _Jpnu, _Tp& _Npnu)
    {
      constexpr auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      constexpr auto _S_inf = __gnu_cxx::__infinity<_Tp>();
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      if (__nu < _Tp{0})
	{
	  _Tp _J_mnu, _N_mnu, _Jp_mnu, _Np_mnu;
	  __cyl_bessel_jn(-__nu, __x, _J_mnu, _N_mnu, _Jp_mnu, _Np_mnu);
	  auto __sinnupi = __sin_pi(-__nu);
	  auto __cosnupi = __cos_pi(-__nu);
	  if (std::abs(__sinnupi) < _S_eps)
	    { // Carefully preserve +-inf.
	      auto __sign = std::copysign(_Tp{1}, __cosnupi);
	      _Jnu = __sign * _J_mnu;
	      _Nnu = __sign * _N_mnu;
	      _Jpnu = __sign * _Jp_mnu;
	      _Npnu = __sign * _Np_mnu;
	    }
	  else if (std::abs(__cosnupi) < _S_eps)
	    { // Carefully preserve +-inf.
	      auto __sign = std::copysign(_Tp{1}, __sinnupi);
	      _Jnu = -__sign * _N_mnu;
	      _Nnu = __sign * _J_mnu;
	      _Jpnu = -__sign * _Np_mnu;
	      _Npnu = __sign * _Jp_mnu;
	    }
	  else
	    {
	      _Jnu = __cosnupi * _J_mnu - __sinnupi * _N_mnu;
	      _Nnu = __sinnupi * _J_mnu + __cosnupi * _N_mnu;
	      _Jpnu = __cosnupi * _Jp_mnu - __sinnupi * _Np_mnu;
	      _Npnu = __sinnupi * _Jp_mnu + __cosnupi * _Np_mnu;
	    }
	}
      else if (__x == _Tp{0})
	{
	  if (__nu == _Tp{0})
	    {
	      _Jnu = _Tp{1};
	      _Jpnu = _Tp{0};
	    }
	  else if (__nu == _Tp{1})
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0.5L};
	    }
	  else
	    {
	      _Jnu = _Tp{0};
	      _Jpnu = _Tp{0};
	    }
	  _Nnu = -_S_inf;
	  _Npnu = _S_inf;
	  return;
	}
      else if (__x > _Tp{1000})
	__cyl_bessel_jn_asymp(__nu, __x, _Jnu, _Nnu, _Jpnu, _Npnu);
      else
	__cyl_bessel_jn_steed(__nu, __x, _Jnu, _Nnu, _Jpnu, _Npnu);
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    void
    __cyl_bessel_jn_neg_arg(_Tp __nu, _Tp __x,
			    std::complex<_Tp>& _Jnu, std::complex<_Tp>& _Nnu,
			    std::complex<_Tp>& _Jpnu, std::complex<_Tp>& _Npnu)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr std::complex<_Tp> _S_i{0, 1};
      if (__x >= _Tp{0})
	std::__throw_domain_error(__N("__cyl_bessel_jn_neg_arg: "
				      "non-negative argument"));
      else
	{
	  _Tp _Jm, _Nm, _Jpm, _Npm;
	  __cyl_bessel_jn(__nu, -__x, _Jm, _Nm, _Jpm, _Npm);
	  auto __phm = std::polar(_Tp{1}, -__nu * _S_pi);
	  auto __php = std::polar(_Tp{1}, __nu * _S_pi);
	  _Jnu = __php * _Jm;
	  _Jpnu = -__php * _Jpm;
	  _Nnu = __phm * _Nm
	       + _S_i * _Tp{2} * __cos_pi(__nu) * _Jm;
	  _Npnu = -__phm * _Npm
	       - _S_i * _Tp{2} * __cos_pi(__nu) * _Jpm;
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
   * @param  __nu  The order of the Bessel function.
   * @param  __x   The argument of the Bessel function.
   * @return  The output Bessel function.
   */
  template<typename _Tp>
    _Tp
    __cyl_bessel_j(_Tp __nu, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__cyl_bessel_j: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__nu >= _Tp{0} && __x * __x < _Tp{10} * (__nu + _Tp{1}))
	return __cyl_bessel_ij_series(__nu, __x, -_Tp{1}, 200);
      else
	{
	  _Tp _J_nu, _N_nu, _Jp_nu, _Np_nu;
	  __cyl_bessel_jn(__nu, __x, _J_nu, _N_nu, _Jp_nu, _Np_nu);
	  return _J_nu;
	}
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
   * @param  __nu  The order of the Neumann function.
   * @param  __x   The argument of the Neumann function.
   * @return  The output Neumann function.
   */
  template<typename _Tp>
    _Tp
    __cyl_neumann_n(_Tp __nu, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__cyl_neumann_n: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else
	{
	  _Tp _J_nu, _N_nu, _Jp_nu, _Np_nu;
	  __cyl_bessel_jn(__nu, __x, _J_nu, _N_nu, _Jp_nu, _Np_nu);
	  return _N_nu;
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
   * @param  __nu  The order of the spherical Neumann function.
   * @param  __x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cyl_hankel_1(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr std::complex<_Tp> _S_i{0, 1};
      if (__nu < _Tp{0})
	return std::exp(-_S_i * _S_pi * __nu) * __cyl_hankel_1(-__nu, __x);
      else if (__isnan(__x))
	return std::complex<_Tp>{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  std::complex<_Tp> _J_n, _N_n, _Jp_n, _Np_n;
	  __cyl_bessel_jn_neg_arg(__nu, __x, _J_n, _N_n, _Jp_n, _Np_n);
	  return _J_n + _S_i * _N_n;
	}
      else
	{
	  _Tp _J_n, _N_n, _Jp_n, _Np_n;
	  __cyl_bessel_jn(__nu, __x, _J_n, _N_n, _Jp_n, _Np_n);
	  return std::complex<_Tp>{_J_n, _N_n};
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
   *   @param  __nu  The order of the spherical Neumann function.
   *   @param  __x  The argument of the spherical Neumann function.
   *   @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __cyl_hankel_2(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      constexpr std::complex<_Tp> _S_i{0, 1};
      if (__nu < _Tp{0})
	return std::exp(_S_i * _S_pi * __nu) * __cyl_hankel_2(-__nu, __x);
      else if (__isnan(__x))
	return std::complex<_Tp>{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  std::complex<_Tp> _J_n, _N_n, _Jp_n, _Np_n;
	  __cyl_bessel_jn_neg_arg(__nu, __x, _J_n, _N_n, _Jp_n, _Np_n);
	  return _J_n - _S_i * _N_n;
	}
      else
	{
	  _Tp _J_n, _N_n, _Jp_n, _Np_n;
	  __cyl_bessel_jn(__nu, __x, _J_n, _N_n, _Jp_n, _Np_n);
	  return std::complex<_Tp>{_J_n, -_N_n};
	}
    }


  /**
   * @brief  Compute the spherical Bessel @f$ j_n(x) @f$
   * 	     and Neumann @f$ n_n(x) @f$ functions and their first
   * 	     derivatives @f$ j_n(x) @f$ and @f$ n'_n(x) @f$
   * 	     respectively.
   *
   * @param  __n  The order of the spherical Bessel function.
   * @param  __x  The argument of the spherical Bessel function.
   * @param[out] __j_n  The output spherical Bessel function.
   * @param[out] __n_n  The output spherical Neumann function.
   * @param[out] __jp_n The output derivative of the spherical Bessel function.
   * @param[out] __np_n The output derivative of the spherical Neumann function.
   */
  template<typename _Tp>
    void
    __sph_bessel_jn(unsigned int __n, _Tp __x,
		    _Tp& __j_n, _Tp& __n_n, _Tp& __jp_n, _Tp& __np_n)
    {
      const auto __nu = _Tp(__n + 0.5L);

      _Tp _J_nu, _N_nu, _Jp_nu, _Np_nu;
      __cyl_bessel_jn(__nu, __x, _J_nu, _N_nu, _Jp_nu, _Np_nu);

      const auto __factor = __gnu_cxx::__math_constants<_Tp>::__root_pi_div_2
			  / std::sqrt(__x);

      __j_n = __factor * _J_nu;
      __n_n = __factor * _N_nu;
      __jp_n = __factor * _Jp_nu - __j_n / (_Tp{2} * __x);
      __np_n = __factor * _Np_nu - __n_n / (_Tp{2} * __x);

      return;
    }

  /**
   * Return the spherical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    void
    __sph_bessel_jn_neg_arg(unsigned int __n, _Tp __x,
			  std::complex<_Tp>& __j_n, std::complex<_Tp>& __n_n,
			  std::complex<_Tp>& __jp_n, std::complex<_Tp>& __np_n)
    {
      if (__x >= _Tp{0})
	std::__throw_domain_error(__N("__sph_bessel_jn_neg_arg: "
				      "non-negative argument"));
      else
	{
	  const auto __nu = _Tp(__n + 0.5L);
	  std::complex<_Tp> _J_nu, _N_nu, _Jp_nu, _Np_nu;
	  __cyl_bessel_jn_neg_arg(__nu, __x, _J_nu, _N_nu, _Jp_nu, _Np_nu);

	  const auto __factor
	    = __gnu_cxx::__math_constants<_Tp>::__root_pi_div_2
	      / std::sqrt(std::complex<_Tp>(__x));

	  __j_n = __factor * _J_nu;
	  __n_n = __factor * _N_nu;
	  __jp_n = __factor * _Jp_nu - __j_n / (_Tp{2} * __x);
	  __np_n = __factor * _Np_nu - __n_n / (_Tp{2} * __x);
	}

      return;
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
   * @param  __n  The non-negative integral order
   * @param  __x  The non-negative real argument
   * @return  The output spherical Bessel function.
   */
  template<typename _Tp>
    _Tp
    __sph_bessel(unsigned int __n, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__sph_bessel: bad argument"));
      else if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__x == _Tp{0})
	{
	  if (__n == 0)
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	{
	  _Tp __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return __j_n;
	}
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
   * @param  __n  The order of the spherical Neumann function.
   * @param  __x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    _Tp
    __sph_neumann(unsigned int __n, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__sph_neumann: bad argument"));
      else if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Tp>();
      else if (__x == _Tp{0})
	return -__gnu_cxx::__infinity<_Tp>();
      else
	{
	  _Tp __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return __n_n;
	}
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
   * @param  __n  The order of the spherical Neumann function.
   * @param  __x  The argument of the spherical Neumann function.
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sph_hankel_1(unsigned int __n, _Tp __x)
    {
      constexpr std::complex<_Tp> _S_i{0, 1};
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      if (__isnan(__x))
	return std::complex<_Tp>{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  std::complex<_Tp> __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn_neg_arg(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return __j_n + _S_i * __n_n;
	}
      else
	{
	  _Tp __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return std::complex<_Tp>{__j_n, __n_n};
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
   * @param  __n  The non-negative integral order
   * @param  __x  The non-negative real argument
   * @return  The output spherical Neumann function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __sph_hankel_2(unsigned int __n, _Tp __x)
    {
      constexpr std::complex<_Tp> _S_i{0, 1};
      constexpr auto _S_nan = __gnu_cxx::__quiet_NaN<_Tp>();
      if (__isnan(__x))
	return std::complex<_Tp>{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  std::complex<_Tp> __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn_neg_arg(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return __j_n - _S_i * __n_n;
	}
      else
	{
	  _Tp __j_n, __n_n, __jp_n, __np_n;
	  __sph_bessel_jn(__n, __x, __j_n, __n_n, __jp_n, __np_n);
	  return std::complex<_Tp>{__j_n, -__n_n};
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_BESSEL_TCC
