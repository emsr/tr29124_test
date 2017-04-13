// Special functions -*- C++ -*-

// Copyright (C) 2006-2017 Free Software Foundation, Inc.
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

/** @file bits/sf_mod_bessel.tcc
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
//     Ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications,
//     Section 9, pp. 355-434, Section 10 pp. 435-478
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 246-249.

#ifndef _GLIBCXX_BITS_SF_MOD_BESSEL_TCC
#define _GLIBCXX_BITS_SF_MOD_BESSEL_TCC 1

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
   * @brief This routine computes the asymptotic modified cylindrical
   * 	    Bessel and functions of order nu: @f$ I_{\nu}(x) @f$,
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
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>
    __cyl_bessel_ik_asymp(_Tp __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__x);
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
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Rsum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
	  _Psum += __ak_xk;
	  const auto __convP = std::abs(__ak_xk) < _S_eps * std::abs(_Psum);

	  ++__k;
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Ssum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
	  _Qsum += __ak_xk;
	  const auto __convQ = std::abs(__ak_xk) < _S_eps * std::abs(_Qsum);

	  if (__convP && __convQ && __k > (__nu / _Tp{2}))
	    break;
	}
      while (__k < _Tp{100} * __nu);

      const auto __coef = std::sqrt(_Tp{1} / (_Tp{2} * _S_pi * __x));
      const auto _Inu = __coef * std::exp(__x) * (_Psum - _Qsum);
      const auto _Ipnu = __coef * std::exp(__x) * (_Rsum - _Ssum);
      const auto _Knu = _S_pi * __coef * std::exp(-__x) * (_Psum + _Qsum);
      const auto _Kpnu =  -_S_pi * __coef * std::exp(-__x) * (_Rsum + _Ssum);

      return __bess_t{__nu, __x, _Inu, _Ipnu, _Knu, _Kpnu};
    }

  /**
   * @brief  Compute the modified Bessel functions @f$ I_\nu(x) @f$ and
   * 	     @f$ K_\nu(x) @f$ and their first derivatives
   * 	     @f$ I'_\nu(x) @f$ and @f$ K'_\nu(x) @f$ respectively.
   * 	     These four functions are computed together for numerical
   * 	     stability.
   *
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>
    __cyl_bessel_ik_steed(_Tp __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_fp_min = _Tp{10} * _S_eps;
      constexpr int _S_max_iter = 15000;
      const auto _S_x_min = _Tp{2};

      const int __n = std::nearbyint(__nu);

      const auto __mu = __nu - _Tp(__n);
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      auto __h = std::max(_S_fp_min, __nu * __xi);
      auto __b = __xi2 * __nu;
      auto __d = _Tp{0};
      auto __c = __h;
      int __i;
      for (__i = 1; __i <= _S_max_iter; ++__i)
	{
	  __b += __xi2;
	  __d = _Tp{1} / (__b + __d);
	  __c = __b + _Tp{1} / __c;
	  const auto __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	}
      if (__i > _S_max_iter)
	return __cyl_bessel_ik_asymp(__nu, __x);

      auto _Inul = _S_fp_min;
      auto _Ipnul = __h * _Inul;
      auto _Inul1 = _Inul;
      auto _Ipnu1 = _Ipnul;
      auto __fact = __nu * __xi;
      for (int __l = __n; __l >= 1; --__l)
	{
	  const auto _Inutemp = __fact * _Inul + _Ipnul;
	  __fact -= __xi;
	  _Ipnul = __fact * _Inutemp + _Inul;
	  _Inul = _Inutemp;
	}

      const auto __f = _Ipnul / _Inul;
      _Tp _Kmu, _Knu1;
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
	  const auto __gamt = __gamma_temme(__mu);
	  auto __ff = __fact
		    * (__gamt.__gamma_1_value * std::cosh(__e)
		     + __gamt.__gamma_2_value * __fact2 * __d);
	  auto __sum = __ff;
	  __e = std::exp(__e);
	  auto __p = __e / (_Tp{2} * __gamt.__gamma_plus_value);
	  auto __q = _Tp{1} / (_Tp{2} * __e * __gamt.__gamma_minus_value);
	  auto __c = _Tp{1};
	  __d = __x2 * __x2;
	  auto __sum1 = __p;
	  int __i;
	  for (__i = 1; __i <= _S_max_iter; ++__i)
	    {
	      __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	      __c *= __d / _Tp(__i);
	      __p /= _Tp(__i) - __mu;
	      __q /= _Tp(__i) + __mu;
	      const auto __del = __c * __ff;
	      __sum += __del;
	      const auto __del1 = __c * (__p - _Tp(__i) * __ff);
	      __sum1 += __del1;
	      if (std::abs(__del) < _S_eps * std::abs(__sum))
		break;
	    }
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__cyl_bessel_ik_steed: "
					   "K-series failed to converge"));
	  _Kmu = __sum;
	  _Knu1 = __sum1 * __xi2;
	}
      else
	{
	  auto __b = _Tp{2} * (_Tp{1} + __x);
	  auto __d = _Tp{1} / __b;
	  auto __delh = __d;
	  auto __h = __delh;
	  auto __q1 = _Tp{0};
	  auto __q2 = _Tp{1};
	  const auto __a1 = _Tp{0.25L} - __mu2;
	  auto __q = __c = __a1;
	  auto __a = -__a1;
	  auto __s = _Tp{1} + __q * __delh;
	  int __i;
	  for (__i = 2; __i <= _S_max_iter; ++__i)
	    {
	      __a -= _Tp{2 * (__i - 1)};
	      __c = -__a * __c / __i;
	      const auto __qnew = (__q1 - __b * __q2) / __a;
	      __q1 = __q2;
	      __q2 = __qnew;
	      __q += __c * __qnew;
	      __b += _Tp{2};
	      __d = _Tp{1} / (__b + __a * __d);
	      __delh = (__b * __d - _Tp{1}) * __delh;
	      __h += __delh;
	      const auto __dels = __q * __delh;
	      __s += __dels;
	      if (std::abs(__dels / __s) < _S_eps)
		break;
	    }
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N("__cyl_bessel_ik_steed: "
					   "Steed's method failed"));
	  __h = __a1 * __h;
	  _Kmu = std::sqrt(_S_pi / (_Tp{2} * __x))
		* std::exp(-__x) / __s;
	  _Knu1 = _Kmu * (__mu + __x + _Tp{0.5L} - __h) * __xi;
	}

      const auto _Kpmu = __mu * __xi * _Kmu - _Knu1;
      const auto _Inumu = __xi / (__f * _Kmu - _Kpmu);
      const auto _Inu = _Inumu * _Inul1 / _Inul;
      const auto _Ipnu = _Inumu * _Ipnu1 / _Inul;
      for (int __i = 1; __i <= __n; ++__i)
	_Kmu = __gnu_cxx::__exchange(_Knu1,
				     (__mu + _Tp(__i)) * __xi2 * _Knu1 + _Kmu);
      const auto _Knu = _Kmu;
      const auto _Kpnu = __nu * __xi * _Kmu - _Knu1;

      return __bess_t{__nu, __x, _Inu, _Ipnu, _Knu, _Kpnu};
    }

  /**
   * @brief  Return the modified cylindrical Bessel functions
   *         and their derivatives of order @f$ \nu @f$ by various means.
   *
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>
    __cyl_bessel_ik(_Tp __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      if (__nu < _Tp{0})
	{
	  const auto _Bessm = __cyl_bessel_ik(-__nu, __x);
	  const auto __arg = -__nu * _S_pi;
	  const auto __sinnupi = __sin_pi(-__nu);
	  if (std::abs(__sinnupi) < _S_eps) // Carefully preserve +-inf.
	    return __bess_t{__nu, __x, _Bessm.__I_value, _Bessm.__I_deriv,
					_Bessm.__K_value, _Bessm.__K_deriv};
	  else
	    return __bess_t{__nu, __x,
	      _Bessm.__I_value + _Tp{2} * __sinnupi * _Bessm.__K_value / _S_pi,
	      _Bessm.__I_deriv + _Tp{2} * __sinnupi * _Bessm.__K_deriv / _S_pi,
	      _Bessm.__K_value, _Bessm.__K_deriv};
	}
      else if (__x == _Tp{0})
	{
	  _Tp _Inu, _Ipnu;
	  if (__nu == _Tp{0})
	    {
	      _Inu = _Tp{1};
	      _Ipnu = _Tp{0};
	    }
	  else if (__nu == _Tp{1})
	    {
	      _Inu = _Tp{0};
	      _Ipnu = _Tp{0.5L};
	    }
	  else
	    {
	      _Inu = _Tp{0};
	      _Ipnu = _Tp{0};
	    }
	  return __bess_t{__nu, __x, _Inu, _Ipnu, _S_inf, -_S_inf};
	}
      else if (__x > _Tp{1000})
	return __cyl_bessel_ik_asymp(__nu, __x);
      else
	return __cyl_bessel_ik_steed(__nu, __x);
    }

  /**
   * @brief  Return the regular modified Bessel function of order
   * 	     @f$ \nu @f$: @f$ I_{\nu}(x) @f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @param  __nu  The order of the regular modified Bessel function.
   * @param  __x   The argument of the regular modified Bessel function.
   * @return  The output regular modified Bessel function.
   */
  template<typename _Tp>
    _Tp
    __cyl_bessel_i(_Tp __nu, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__cyl_bessel_i: bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__nu >= _Tp{0} && __x * __x < _Tp{10} * (__nu + _Tp{1}))
	return __cyl_bessel_ij_series(__nu, __x, +_Tp{1}, 200);
      else
	return __cyl_bessel_ik(__nu, __x).__I_value;
    }

  /**
   * @brief  Return the irregular modified Bessel function
   * 	     @f$ K_{\nu}(x) @f$ of order @f$ \nu @f$.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @param  __nu  The order of the irregular modified Bessel function.
   * @param  __x   The argument of the irregular modified Bessel function.
   * @return  The output irregular modified Bessel function.
   */
  template<typename _Tp>
    _Tp
    __cyl_bessel_k(_Tp __nu, _Tp __x)
    {
      if (__x < _Tp{0})
	std::__throw_domain_error(__N("__cyl_bessel_k: Bad argument"));
      else if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else
	return __cyl_bessel_ik(__nu, __x).__K_value;
    }

  /**
   * @brief  Compute the spherical modified Bessel functions
   * 	     @f$ i_n(x) @f$ and @f$ k_n(x) @f$ and their first
   * 	     derivatives @f$ i'_n(x) @f$ and @f$ k'_n(x) @f$
   * 	     respectively.
   *
   * @param  __n  The order of the modified spherical Bessel function.
   * @param  __x  The argument of the modified spherical Bessel function.
   * @return A struct containing the modified spherical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__sph_mod_bessel_t<unsigned int, _Tp, _Tp>
    __sph_bessel_ik(unsigned int __n, _Tp __x)
    {
      using __sph_t = __gnu_cxx::__sph_mod_bessel_t<unsigned int, _Tp, _Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);

      if (__isnan(__x))
	return __sph_t{__n, __x, _S_NaN, _S_NaN, _S_NaN, _S_NaN};
      else if (__x == _Tp{0})
	{
	  const auto _S_inf = __gnu_cxx::__infinity(__x);
	  if (__n == 0)
	    return __sph_t{__n, __x, _Tp{1}, _Tp{0}, _S_inf, -_S_inf};
	  else
	    return __sph_t{__n, __x, _Tp{0}, _Tp{0}, _S_inf, -_S_inf};
	}
      else
	{
	  const auto __nu = _Tp(__n + 0.5L);
	  auto _Bess = __cyl_bessel_ik(__nu, __x);

	  const auto __factor = __gnu_cxx::__const_root_pi_div_2(_Tp{})
			      / std::sqrt(__x);

	  const auto __i_n = __factor * _Bess.__I_value;
	  const auto __ip_n = __factor * _Bess.__I_deriv
			    - __i_n / (_Tp{2} * __x);
	  const auto __k_n = __factor * _Bess.__K_value;
	  const auto __kp_n = __factor * _Bess.__K_deriv
			    - __k_n / (_Tp{2} * __x);

	  return __sph_t{__n, __x, __i_n, __ip_n, __k_n, __kp_n};
	}
    }


  /**
   * @brief  Compute the Airy functions
   * 	     @f$ Ai(x) @f$ and @f$ Bi(x) @f$ and their first
   * 	     derivatives @f$ Ai'(x) @f$ and @f$ Bi(x) @f$
   * 	     respectively.
   *
   * @param  __z  The argument of the Airy functions.
   * @return A struct containing the Airy functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__airy_t<_Tp, _Tp>
    __airy(_Tp __z)
    {
      using __ai_t = __gnu_cxx::__airy_t<_Tp, _Tp>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__z);
      const auto _S_inf = __gnu_cxx::__infinity(__z);
      const auto _S_pi = __gnu_cxx::__const_pi(__z);
      const auto _S_sqrt3 = __gnu_cxx::__const_root_3(__z);
      const auto __absz = std::abs(__z);
      const auto __rootz = std::sqrt(__absz);
      const auto __xi = _Tp{2} * __absz * __rootz / _Tp{3};

      if (__isnan(__z))
	return __ai_t{__z, _S_NaN, _S_NaN, _S_NaN, _S_NaN};
      else if (__z == _S_inf)
	return __ai_t{__z, _Tp{0}, _Tp{0}, _S_inf, _S_inf};
      else if (__z == -_S_inf)
	return __ai_t{__z, _Tp{0}, _Tp{0}, _Tp{0}, _Tp{0}};
      else if (__z > _Tp{0})
	{
	  const auto _Bess13 = __cyl_bessel_ik(_Tp{1} / _Tp{3}, __xi);
	  const auto _Ai = __rootz * _Bess13.__K_value / (_S_sqrt3 * _S_pi);
	  const auto _Bi = __rootz * (_Bess13.__K_value / _S_pi
				    + _Tp{2} * _Bess13.__I_value / _S_sqrt3);

	  const auto _Bess23 = __cyl_bessel_ik(_Tp{2} / _Tp{3}, __xi);
	  const auto _Aip = -__z * _Bess23.__K_value / (_S_sqrt3 * _S_pi);
	  const auto _Bip = __z * (_Bess23.__K_value / _S_pi
				 + _Tp{2} * _Bess23.__I_value / _S_sqrt3);

	  return __ai_t{__z, _Ai, _Aip, _Bi, _Bip};
	}
      else if (__z < _Tp{0})
	{
	  const auto _Bess13 = __cyl_bessel_jn(_Tp{1} / _Tp{3}, __xi);
	  const auto _Ai = +__rootz * (_Bess13.__J_value
				     - _Bess13.__N_value / _S_sqrt3) / _Tp{2};
	  const auto _Bi = -__rootz * (_Bess13.__N_value
				     + _Bess13.__J_value / _S_sqrt3) / _Tp{2};

	  const auto _Bess23 = __cyl_bessel_jn(_Tp{2} / _Tp{3}, __xi);
	  const auto _Aip = __absz * (_Bess23.__N_value / _S_sqrt3
				    + _Bess23.__J_value) / _Tp{2};
	  const auto _Bip = __absz * (_Bess23.__J_value / _S_sqrt3
				    - _Bess23.__N_value) / _Tp{2};

	  return __ai_t{__z, _Ai, _Aip, _Bi, _Bip};
	}
      else
	{
	  // Reference:
	  //  Abramowitz & Stegun, page 446 section 10.4.4 on Airy functions.
	  // The number is Ai(0) = 3^{-2/3}/\Gamma(2/3).
	  const auto _Ai
	    = _Tp{0.3550280538878172392600631860041831763979791741991772L};
	  const auto _Bi = _Ai * _S_sqrt3;

	  // Reference:
	  //  Abramowitz & Stegun, page 446 section 10.4.5 on Airy functions.
	  // The number is Ai'(0) = -3^{-1/3}/\Gamma(1/3).
	  const auto _Aip
	    = -_Tp{0.25881940379280679840518356018920396347909113835493L};
	  const auto _Bip = -_Aip * _S_sqrt3;

	  return __ai_t{__z, _Ai, _Aip, _Bi, _Bip};
	}
    }

  /**
   * @brief  Compute the Fock-type Airy functions
   * 	     @f$ w_1(x) @f$ and @f$ w_2(x) @f$ and their first
   * 	     derivatives @f$ w_1'(x) @f$ and @f$ w_2'(x) @f$
   * 	     respectively.
   * @f[
   *   w_1(x) = \sqrt{\pi}(Ai(x) + iBi(x))
   * @f]
   * @f[
   *   w_2(x) = \sqrt{\pi}(Ai(x) - iBi(x))
   * @f]
   *
   * @param  __x   The argument of the Airy functions.
   * @return A struct containing the Fock-type Airy functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__fock_airy_t<_Tp, std::complex<_Tp>>
    __fock_airy(_Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      using __fock_t = __gnu_cxx::__fock_airy_t<_Tp, _Cmplx>;
      const auto _S_sqrtpi = __gnu_cxx::__const_root_pi(__x);

      const auto _Ai = __airy(__x);

      const auto __w1 = _S_sqrtpi * _Cmplx(_Ai.__Ai_value, _Ai.__Bi_value);
      const auto __w1p = _S_sqrtpi * _Cmplx(_Ai.__Ai_deriv, _Ai.__Bi_deriv);
      const auto __w2 = _S_sqrtpi * _Cmplx(_Ai.__Ai_value, -_Ai.__Bi_value);
      const auto __w2p = _S_sqrtpi * _Cmplx(_Ai.__Ai_deriv, -_Ai.__Bi_deriv);

      return __fock_t{__x, __w1, __w1p, __w2, __w2p};
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_MOD_BESSEL_TCC
