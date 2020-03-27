// Special functions -*- C++ -*-

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

#ifndef _GLIBCXX_BITS_SF_BESSEL_TCC
#define _GLIBCXX_BITS_SF_BESSEL_TCC 1

#pragma GCC system_header

#include <ext/continued_fractions.h>

#define BESSEL_DEBUG 1
#ifdef BESSEL_DEBUG
#  include <iostream>
#  include <iomanip>
#endif

namespace std _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

namespace __detail
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
   * @param  __nu  The order of the Bessel function.
   * @param  __x   The argument of the Bessel function.
   * @param  __sgn  The sign of the alternate terms
   * 		    -1 for the Bessel function of the first kind.
   * 		    +1 for the modified Bessel function of the first kind.
   * @param  __max_iter  The maximum number of iterations for sum.
   * @return  The output Bessel function.
   */
  template<typename _Tnu, typename _Tp>
    constexpr _Tp
    __cyl_bessel_ij_series(_Tnu __nu, _Tp __x, int __sgn,
			   unsigned int __max_iter)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      if (std::abs(__x) < _S_eps)
	{
	  if (__nu == _Tnu{0})
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	{
	  const auto __x2 = __x / _Real{2};

	  const auto __xx4 = _Val(__sgn) * __x2 * __x2;
	  auto _Jn = _Val{1};
	  auto __term = _Val{1};
	  for (unsigned int __i = 1; __i < __max_iter; ++__i)
	    {
	      __term *= __xx4 / (_Val(__i) * (_Val(__nu) + _Val(__i)));
	      _Jn += __term;
	      if (std::abs(__term / _Jn) < _S_eps)
		break;
	    }

	  auto __fact = _Val(__nu) * std::log(__x2);
	  __fact -= __log_gamma(_Real{1} + __nu);
	  __fact = std::exp(__fact);
	  return __fact * _Jn;
	}
    }

  /**
   * A type for Bessel asymptotic sums.
   */
  template<typename _Tnu, typename _Tp>
    struct __cyl_bessel_asymp_sums_t
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
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
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename _Tp>
    constexpr __cyl_bessel_asymp_sums_t<_Tnu, _Tp>
    __cyl_bessel_asymp_sums(_Tnu __nu, _Tp __x, int __sgn)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using __bess_t = __cyl_bessel_asymp_sums_t<_Tnu, _Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon<_Real>();
      const auto __2nu = _Real{2} * __nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __r8x = _Tp{1} / (_Real{8} * __x);
      const auto __nu_min = std::real(__nu / _Real{2});
      const auto __nu_max = std::abs(_Real{100} * (__nu + _Tnu{1}));
      auto __k = 0;
      auto __bk_xk = _Val{1};
      auto _Rsum = __bk_xk;
      auto __ak_xk = _Val{1};
      auto _Psum = __ak_xk;
      auto __convP = false;
      ++__k;
      auto __2km1 = 1;
      __bk_xk *= (__4nu2 + 3) * __r8x;
      auto _Ssum = __bk_xk;
      __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * __r8x;
      auto _Qsum = __ak_xk;
      auto __convQ = false;
      auto __ak_xk_prev = std::abs(__ak_xk);
      do
	{
	  ++__k;
	  auto __rk8x = __r8x / _Real(__k);
	  __2km1 += 2;
	  __bk_xk = __sgn * (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk * __rk8x;
	  _Rsum += __bk_xk;
	  __ak_xk *= __sgn * (__2nu - __2km1) * (__2nu + __2km1) * __rk8x;
	  if (__k > __nu_min && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Psum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convP = std::abs(__ak_xk) < _S_eps * std::abs(_Psum);

	  ++__k;
	  __rk8x = __r8x / _Real(__k);
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk * __rk8x;
	  _Ssum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) * __rk8x;
	  if (__k > __nu_min && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Qsum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convQ = std::abs(__ak_xk) < _S_eps * std::abs(_Qsum);

	  if (__convP && __convQ)
	    break;
	}
      while (__k < __nu_max);

      return __bess_t{_Psum, _Qsum, _Rsum, _Ssum};
    }

  /**
   *
   */
  template<typename _Tnu, typename _Tp>
    constexpr __gnu_cxx::__cyl_bessel_t<_Tnu, _Tp, _Tp>
    __cyl_bessel_jn_asymp(_Tnu __nu, _Tp __x)
    {
      // FIXME: This will promote float to double if _Tnu is integral.
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using __bess_t = __gnu_cxx::__cyl_bessel_t<_Tnu, _Tp, _Tp>;
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Real>;
      const auto _S_pi_2 = __gnu_cxx::numbers::__pi_half_v<_Real>;

      const auto __sums = __cyl_bessel_asymp_sums(__nu, __x, -1);

      const auto __omega = __x - (__nu + _Real{0.5L}) * _S_pi_2;
      const auto __c = std::cos(__omega);
      const auto __s = std::sin(__omega);

      const auto __coef = std::sqrt(_Real{2} / (_S_pi * __x));
      return __bess_t{__nu, __x,
		 __coef * (__c * __sums._Psum - __s * __sums._Qsum),
		-__coef * (__s * __sums._Rsum + __c * __sums._Ssum),
		 __coef * (__s * __sums._Psum + __c * __sums._Qsum),
		 __coef * (__c * __sums._Rsum - __s * __sums._Ssum)};
    }

  /**
   * Compute ratios of Bessel functions using the S-fraction.
   *
   * @param  __nu  The order of the Hankel ratio.
   * @param  __x   The argument of the Hankel ratio.
   * @param  __zeta A variable encapsulating the regular and irregular
   *           Bessel functions; Use @f$ \zeta = (iz)^2 @f$ for @f$ J_\nu(z) @f$
   *           and @f$ \zeta = z^2 @f$ for @f$ I_\nu(z) @f$.
   */
  template<typename _Tnu, typename _Tp, typename _Tzeta>
    std::complex<__gnu_cxx::__num_traits_t<
		 __gnu_cxx::fp_promote_t<_Tnu, _Tp, _Tzeta>>>
    __cyl_bessel_ratio_s_frac(_Tnu __nu, _Tp __z, _Tzeta __zeta)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      auto __a_J
	= [__nu, __z, __zeta](std::size_t __k, _Tp)
	  {
	    using __type = decltype(_Tnu{} * _Tzeta{});
	    if (__k == 1)
	      return __type(__z / (_Tnu{2} * __nu + _Tnu{2}));
	    else
	      return __zeta
		   / (_Tnu{4} * (__nu + _Tnu(__k - 1)) * (__nu + _Tnu(__k)));
	  };
      using _AFun = decltype(__a_J);

      auto __b_J = [](std::size_t, _Tp) -> _Real { return _Real{1}; };
      using _BFun = decltype(__b_J);

      auto __w_J = [](std::size_t, _Tp) -> _Real { return _Real{0}; };
      using _WFun = decltype(__w_J);

      _SteedContinuedFraction<_Tp, _AFun, _BFun, _WFun>
      _J(__a_J, __b_J, __w_J);

      // b_0 is 0 not 1 so subtract 1.
      return _J(__z) - _Real{1};
    }

  /**
   * 
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<__gnu_cxx::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_j_ratio_s_frac(_Tnu __nu, _Tp __z)
    {
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __iz = _Cmplx{0, 1} * __z;
      const auto __zeta = __iz * __iz;
      const auto _Jrat = __cyl_bessel_ratio_s_frac(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Jrat);
      else
	return _Jrat;
    }

  /**
   * 
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<__gnu_cxx::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_i_ratio_s_frac(_Tnu __nu, _Tp __z)
    {
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;

      const auto __zeta = __z * __z;
      const auto _Irat = __cyl_bessel_ratio_s_frac(__nu, __z, __zeta);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
	return std::real(_Irat);
      else
	return _Irat;
    }

  /**
   * Compute ratios of Hankel functions using the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    std::complex<__gnu_cxx::__num_traits_t<
		 __gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_ratio_j_frac(_Tnu __nu, _Tp __z, _Tp __sgn)
    {
      using _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>;
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto __zeta = _Cmplx{0, 2} * __z;
      using _Tzeta = decltype(__zeta);

      auto __a_H
	= [__nu](std::size_t __k, _Tp)
	  {
	    const auto __kk = _Tnu(2 * __k - 1) / _Tnu{2};
	    return (__nu - __kk) * (__nu + __kk);
	  };
      using _NumFun = decltype(__a_H);

      auto __b_H
	= [__zeta, __sgn](std::size_t __k, _Tp)
	  { return __sgn * _Tzeta(2 * __k) + __zeta; };
      using _DenFun = decltype(__b_H);

      auto __w_H
	= [__zeta](std::size_t __k, _Tp)
	  { return _Tzeta(__k) + __zeta / _Tzeta{2}; };
      using _TailFun = decltype(__w_H);

      _SteedContinuedFraction<_Tp, _NumFun, _DenFun, _TailFun>
      _H(__a_H, __b_H, __w_H);

      return (_Tzeta(2 * __nu + 1) + __sgn * __zeta) / (_Tp{2} * __z)
	   + __sgn * (_H(__z) - __b_H(0, _Tp{})) / __z;
    }

  /**
   * Return the Hankel function ratio of the first kind from the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    inline std::complex<__gnu_cxx::__num_traits_t<
			__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_1_ratio_j_frac(_Tnu __nu, _Tp __z)
    { return __cyl_hankel_ratio_j_frac(__nu, __z, _Tp{-1}); }

  /**
   * Return the Hankel function ratio of the second kind from the J-fraction.
   */
  template<typename _Tnu, typename _Tp>
    inline std::complex<__gnu_cxx::__num_traits_t<
			__gnu_cxx::fp_promote_t<_Tnu, _Tp>>>
    __cyl_hankel_2_ratio_j_frac(_Tnu __nu, _Tp __z)
    { return __cyl_hankel_ratio_j_frac(__nu, __z, _Tp{+1}); }

  /**
   * Return the modified Bessel function ratio of the second kind
   * from the J-fraction ratios of Hankel functions.
   */
  template<typename _Tnu, typename _Tp,
	   typename _Val = __gnu_cxx::fp_promote_t<_Tnu, _Tp>>
    std::conditional_t<__gnu_cxx::is_complex_v<_Val>,
			std::complex<__gnu_cxx::__num_traits_t<_Val>>,
			_Val>
    __cyl_bessel_k_ratio_j_frac(_Tnu __nu, _Tp __z)
    {
      using _Real = __gnu_cxx::__num_traits_t<_Val>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Real>;
      const auto __ph = std::arg(__z);

      _Cmplx _Krat;
      if (__ph > -_S_pi && __ph <= _S_pi / _Real{2})
	_Krat = _S_i * __cyl_hankel_1_ratio_j_frac(__nu, _S_i * __z);
      else
	_Krat = -_S_i * __cyl_hankel_2_ratio_j_frac(__nu, -_S_i * __z);

      if constexpr (!__gnu_cxx::is_complex_v<_Val>)
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
   * @param __mu The input parameter of the gamma functions.
   * @return  An output structure containing four gamma functions.
   */
  template<typename _Tp>
    __gnu_cxx::__gamma_temme_t<_Tp>
    __gamma_temme(_Tp __mu)
    {
      using __gammat_t = __gnu_cxx::__gamma_temme_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__mu);
      const auto _S_gamma_E = __gnu_cxx::numbers::__gamma_e_v<_Tp>;

      if (std::abs(__mu) < _S_eps)
	return __gammat_t{__mu, _Tp{1}, _Tp{1}, -_S_gamma_E, _Tp{1}};
      else
	{
	  _Tp __gamp, __gamm;
	  if (std::real(__mu) <= _Tp{0})
	    {
	      __gamp = __gamma_reciprocal_series(_Tp{1} + __mu);
	      __gamm = -__gamma_reciprocal_series(-__mu) / __mu;
	    }
	  else
	    {
	      __gamp = __gamma_reciprocal_series(__mu) / __mu;
	      __gamm = __gamma_reciprocal_series(_Tp{1} - __mu);
	    }
	  const auto __gam1 = (__gamm - __gamp) / (_Tp{2} * __mu);
	  const auto __gam2 = (__gamm + __gamp) / _Tp{2};
	  return __gammat_t{__mu, __gamp, __gamm, __gam1, __gam2};
	}
    }

  /**
   * A little return type for the Temme series.
   */
  template<typename _Tp>
    struct __bessel_nk_series_t
    {
      _Tp _Z_mu;
      _Tp _Z_mup1;
    };

  /**
   * This routine computes the dominant cylindrical bessel function solutions
   * by series summation for order @f$ |\mu| < 1/2 @f$.
   *
   * @param __mu The order of the Bessel functions @f$ |\mu| < 1/2 @f$.
   * @param __x  The argument of the Bessel functions.
   * @param __modified If true solve for the modified Bessel function
   *                   @f$ K_\mu @f$, otherwise solve for the Neumann function
   *                   @f$ N_\mu @f$,
   * @return A structure containing Z_{\mu} and Z_{\mu+1}.
   */
  template<typename _Tp>
    __bessel_nk_series_t<_Tp>
    __cyl_bessel_nk_series(_Tp __mu, _Tp __x, bool __modified = false,
			   int __max_iter = 100)
    {
      const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      const auto __xi = _Tp{1} / __x;
      const auto __x2 = __x / _Tp{2};

      const auto __fact = _Tp{1} / std::__detail::__sinc_pi(__mu);
      const auto __lx2 = -std::log(__x2);
      const auto __arg = __mu * __lx2;
      const auto __fact2 = std::__detail::__sinhc(__arg);
      const auto __gamt = std::__detail::__gamma_temme(__mu);
      const auto __norm = __modified ? _Tp{-1} : _Tp{2} / _S_pi;
      auto __ff = __norm * __fact
		* (__gamt.__gamma_1_value * std::cosh(__arg)
		 + __gamt.__gamma_2_value * __fact2 * __lx2);
      const auto __e = std::exp(__arg);
      auto __p = __norm * __e / (_Tp{2} * __gamt.__gamma_plus_value);
      auto __q = __norm / (__e * _Tp{2} * __gamt.__gamma_minus_value);
      const auto __fact3 = __modified
			 ? _Tp{0}
			 : std::__detail::__sinc_pi(__mu / _Tp{2});
      const auto __r = __modified
		     ? _Tp{0}
		     : __fact3 * __fact3 * _S_pi * _S_pi * __mu / _Tp{2};
      auto __c = _Tp{1};
      const auto __d = __modified ? __x2 * __x2 : -__x2 * __x2;
      auto __sum_mu = __ff + __r * __q;
      auto __sum_mup1 = __p;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __ff = (__i * __ff + __p + __q)
	       / ((_Tp(__i) - __mu) * (_Tp(__i) + __mu));
	  __c *= __d / _Tp(__i);
	  __p /= _Tp(__i) - __mu;
	  __q /= _Tp(__i) + __mu;
	  const auto __del_mu = __c * (__ff + __r * __q);
	  __sum_mu += __del_mu;
	  const auto __del_mup1 = __c * __p - _Tp(__i) * __del_mu;
	  __sum_mup1 += __del_mup1;
	  if (std::abs(__del_mu) < _S_eps * std::abs(__sum_mu))
	    break;
	}
      if (__i > __max_iter)
	std::__throw_runtime_error(__N("__cyl_bessel_nk_series: "
				       "Series failed to converge"));
      auto _N_mu = -__sum_mu;
      auto _N_mup1 = -_Tp{2} * __xi * __sum_mup1;

      return {_N_mu, _N_mup1};
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
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp>
    __cyl_bessel_jn_steed(_Tp __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_tiny = __gnu_cxx::__lim_min(__x);
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      // When the multiplier is N i.e.
      // fp_min = N * min()
      // Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
      //const _Tp _S_fp_min = _Tp{20} * __gnu_cxx::__lim_min(__nu);
      constexpr int _S_max_iter = 15000;
      const auto _S_x_min = _Tp{2};
      const auto _S_fp_min = __gnu_cxx::__sqrt_min(__nu);

      const int __n = (__x < _S_x_min
		    ? std::nearbyint(__nu)
		    : std::max(0,
			       static_cast<int>(__nu - __x + _Tp{1.5L})));

      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      const auto _Wronski = __xi2 / _S_pi;
      int __isign = 1;
//#if 1
      auto __h = std::max(_S_fp_min, __nu * __xi);
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
	  __d = _Tp{1} / __d;
	  __c = __b - _Tp{1} / __c;
	  if (std::abs(__c) < _S_fp_min)
	    __c = _S_fp_min;
	  const auto __del = __c * __d;
	  __h *= __del;
	  if (__d < _Tp{0})
	    __isign = -__isign;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	}
      if (__i > _S_max_iter)
	return __cyl_bessel_jn_asymp(__nu, __x);
//#else
      _Tp __h2;
      try
	{
	  __h2 = __nu * __xi - __cyl_bessel_j_ratio_s_frac(__nu, __x);
	}
      catch (...)
	{
	  return __cyl_bessel_jn_asymp(__nu, __x);
	}
//#endif
#ifdef BESSEL_DEBUG
std::cerr << ' ' << std::setw(10) << __nu
	  << ' ' << std::setw(10) << __x
	  << ' ' << std::setw(2) << __isign
	  << ' ' << std::setw(10) << __h2
	  << ' ' << std::setw(10) << __h
	  << ' ' << std::setw(14) << __h2 - __h
	  << '\n';
#endif

      // Recur recessive solution downward to |mu| < 1/2.
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

      const auto __mu = __nu - _Tp(__n);
      const auto __mu2 = __mu * __mu;
      const auto __f = _Jpnul / _Jnul;
      _Tp _Nmu, _Nnu1, _Npmu, _Jmu;
      if (__x < _S_x_min)
	{
	  const auto _Z = __cyl_bessel_nk_series(__mu, __x);
	  _Nmu = _Z._Z_mu;
	  _Nnu1 = _Z._Z_mup1;
	  _Npmu = __mu * __xi * _Nmu - _Nnu1;
	  _Jmu = _Wronski / (_Npmu - __f * _Nmu);
	}
      else
	{
	  auto __pq = __cyl_hankel_1_ratio_j_frac(__mu, __x);
	  auto [__p, __q] = reinterpret_cast<_Tp(&)[2]>(__pq);
	  // FIXME: Is this right?
	  __q = std::abs(__q);
	  const auto __gam = (__p - __f) / __q;
	  _Jmu = std::sqrt(_Wronski / ((__p - __f) * __gam + __q));
	  _Jmu = std::copysign(_Jmu, _Jnul);
	  _Nmu = __gam * _Jmu;
	  _Npmu = (__p + __q / __gam) * _Nmu;
	  _Nnu1 = __mu * __xi * _Nmu - _Npmu;
        }

      __fact = _Jmu / _Jnul;
      const auto _Jnu = __fact * _Jnul1;
      const auto _Jpnu = __fact * _Jpnu1;
      if (std::abs(_S_pi * __x * _Jnu / _Tp{2}) > _S_tiny)
	{
	  for (int __i = 1; __i <= __n; ++__i)
	    _Nmu = __gnu_cxx::__exchange(_Nnu1,
					 (__mu + __i) * __xi2 * _Nnu1 - _Nmu);
	  const auto _Nnu = _Nmu;
	  const auto _Npnu = __nu * __xi * _Nmu - _Nnu1;
	  return __bess_t{__nu, __x, _Jnu, _Jpnu, _Nnu, _Npnu};
	}
      else
	return __bess_t{__nu, __x, _Jnu, _Jpnu, -_S_inf, _S_inf};
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   * of order @f$ \nu @f$ by various means.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp>
    __cyl_bessel_jn(_Tp __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_inf = __gnu_cxx::__infinity(__x);
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      if (__nu < _Tp{0})
	{
	  const auto _Bess = __cyl_bessel_jn(-__nu, __x);
	  const auto __sinnupi = __sin_pi(-__nu);
	  const auto __cosnupi = __cos_pi(-__nu);
	  if (std::abs(__sinnupi) < _S_eps)
	    { // Carefully preserve +-inf.
	      const auto __sign = std::copysign(_Tp{1}, __cosnupi);
	      return __bess_t{__nu, __x,
			__sign * _Bess.__J_value, __sign * _Bess.__J_deriv,
			__sign * _Bess.__N_value, __sign * _Bess.__N_deriv};
	    }
	  else if (std::abs(__cosnupi) < _S_eps)
	    { // Carefully preserve +-inf.
	      const auto __sign = std::copysign(_Tp{1}, __sinnupi);
	      return __bess_t{__nu, __x,
			-__sign * _Bess.__N_value, -__sign * _Bess.__N_deriv,
			 __sign * _Bess.__J_value,  __sign * _Bess.__J_deriv};
	    }
	  else
	    {
	      return __bess_t{__nu, __x,
		__cosnupi * _Bess.__J_value - __sinnupi * _Bess.__N_value,
		__cosnupi * _Bess.__J_deriv - __sinnupi * _Bess.__N_deriv,
		__sinnupi * _Bess.__J_value + __cosnupi * _Bess.__N_value,
		__sinnupi * _Bess.__J_deriv + __cosnupi * _Bess.__N_deriv};
	    }
	}
      else if (__x == _Tp{0})
	{
	  _Tp _Jnu, _Jpnu;
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
	  return __bess_t{__nu, __x, _Jnu, _Jpnu, -_S_inf, _S_inf};
	}
      else if (__x > _Tp{1000})
	return __cyl_bessel_jn_asymp(__nu, __x);
      else
	return __cyl_bessel_jn_steed(__nu, __x);
    }

  /**
   * @brief  Return the cylindrical Bessel functions and their derivatives
   *         of real order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, std::complex<_Tp>>
    __cyl_bessel_jn_neg_arg(_Tp __nu, _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      using __bess_t = __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Cmplx>;
      constexpr _Cmplx _S_i{0, 1};
      if (__x >= _Tp{0})
	std::__throw_domain_error(__N("__cyl_bessel_jn_neg_arg: "
				      "non-negative argument"));
      else
	{
	  const auto _Bess = __cyl_bessel_jn(__nu, -__x);
	  const auto __phm = __polar_pi(_Tp{1}, -__nu);
	  const auto __php = __polar_pi(_Tp{1}, __nu);
	  const auto __cosp = __cos_pi(__nu);
	  return __bess_t{__nu, __x,
			  __php * _Bess.__J_value,
			  -__php * _Bess.__J_deriv,
			  __phm * _Bess.__N_value
				+ _S_i * _Tp{2} * __cosp * _Bess.__J_value,
			  -__phm * _Bess.__N_deriv
				- _S_i * _Tp{2} * __cosp * _Bess.__J_deriv};
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
      else if (std::isnan(__nu) || std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__nu >= _Tp{0} && __x * __x < _Tp{10} * (__nu + _Tp{1}))
	return __cyl_bessel_ij_series(__nu, __x, -1, 200);
      else
	return __cyl_bessel_jn(__nu, __x).__J_value;
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
      else if (std::isnan(__nu) || std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else
	return __cyl_bessel_jn(__nu, __x).__N_value;
    }

  /**
   * @brief  Return the cylindrical Hankel functions of the first
   *         and second kinds and their derivatives.
   *
   */
  template<typename _Tp>
    __gnu_cxx::__cyl_hankel_t<_Tp, _Tp, std::complex<_Tp>>
    __cyl_hankel_h1h2(_Tp __nu, _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx _S_i{0, 1};

      _Cmplx __ph1 = _Tp{1}, __ph2 = _Tp{1};
      if (__nu < _Tp{0})
	{
	  __ph1 = __polar_pi(_Tp{1}, -__nu);
	  __ph2 = __polar_pi(_Tp{1}, +__nu);
	  __nu = -__nu;
	}

      // The two _Bess types are different.
      // We might still be able to assign the real output to the complex one.
      if (__x < _Tp{0})
	{
	  const auto _Bess = __cyl_bessel_jn_neg_arg(__nu, __x);
	  const auto _H1 = __ph1 * (_Bess.__J_value + _S_i * _Bess.__N_value);
	  const auto _dH1 = __ph1 * (_Bess.__J_deriv + _S_i * _Bess.__N_deriv);
	  const auto _H2 = __ph2 * (_Bess.__J_value - _S_i * _Bess.__N_value);
	  const auto _dH2 = __ph2 * (_Bess.__J_deriv - _S_i * _Bess.__N_deriv);
	  return {__nu, __x, _H1, _dH1, _H2, _dH2};
	}
      else
	{
	  const auto _Bess = __cyl_bessel_jn(__nu, __x);
	  const auto _H1 = __ph1 * _Cmplx{_Bess.__J_value, _Bess.__N_value};
	  const auto _dH1 = __ph1 * _Cmplx{_Bess.__J_deriv, _Bess.__N_deriv};
	  const auto _H2 = __ph2 * _Cmplx{_Bess.__J_value, -_Bess.__N_value};
	  const auto _dH2 = __ph2 * _Cmplx{_Bess.__J_deriv, -_Bess.__N_deriv};
	  return {__nu, __x, _H1, _dH1, _H2, _dH2};
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
      using _Cmplx = std::complex<_Tp>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__x);
      constexpr _Cmplx _S_i{0, 1};
      if (__nu < _Tp{0})
	return __polar_pi(_Tp{1}, -__nu)
	     * __cyl_hankel_1(-__nu, __x);
      else if (std::isnan(__x))
	return _Cmplx{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  const auto _Bess = __cyl_bessel_jn_neg_arg(__nu, __x);
	  return _Bess.__J_value + _S_i * _Bess.__N_value;
	}
      else
	{
	  const auto _Bess = __cyl_bessel_jn(__nu, __x);
	  return _Cmplx{_Bess.__J_value, _Bess.__N_value};
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
      using _Cmplx = std::complex<_Tp>;
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__x);
      constexpr _Cmplx _S_i{0, 1};
      if (__nu < _Tp{0})
	return __polar_pi(_Tp{1}, __nu)
	     * __cyl_hankel_2(-__nu, __x);
      else if (std::isnan(__x))
	return _Cmplx{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  const auto _Bess = __cyl_bessel_jn_neg_arg(__nu, __x);
	  return _Bess.__J_value - _S_i * _Bess.__N_value;
	}
      else
	{
	  const auto _Bess = __cyl_bessel_jn(__nu, __x);
	  return _Cmplx{_Bess.__J_value, -_Bess.__N_value};
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
   * @return  The output derivative of the spherical Neumann function.
   */
  template<typename _Tp>
    __gnu_cxx::__sph_bessel_t<unsigned int, _Tp, _Tp>
    __sph_bessel_jn(unsigned int __n, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__sph_bessel_t<unsigned int, _Tp, _Tp>;
      const auto __nu = _Tp(__n + 0.5L);

      const auto _Bess = __cyl_bessel_jn(__nu, __x);

      const auto __factor = __gnu_cxx::numbers::__root_pi_div_2_v<_Tp>
			  / std::sqrt(__x);

      const auto __j_n = __factor * _Bess.__J_value;
      const auto __jp_n = __factor * _Bess.__J_deriv - __j_n / (_Tp{2} * __x);
      const auto __n_n = __factor * _Bess.__N_value;
      const auto __np_n = __factor * _Bess.__N_deriv - __n_n / (_Tp{2} * __x);

      return __bess_t{__n, __x, __j_n, __jp_n, __n_n, __np_n};
    }

  /**
   * Return the spherical Bessel functions and their derivatives
   * of order @f$ \nu @f$ and argument @f$ x < 0 @f$.
   */
  template<typename _Tp>
    __gnu_cxx::__sph_bessel_t<unsigned int, _Tp, std::complex<_Tp>>
    __sph_bessel_jn_neg_arg(unsigned int __n, _Tp __x)
    {
      using _Cmplx = std::complex<_Tp>;
      using __bess_t
	= __gnu_cxx::__sph_bessel_t<unsigned int, _Tp, _Cmplx>;
      if (__x >= _Tp{0})
	std::__throw_domain_error(__N("__sph_bessel_jn_neg_arg: "
				      "non-negative argument"));
      else
	{
	  const auto __nu = _Tp(__n + 0.5L);
	  const auto _Bess = __cyl_bessel_jn_neg_arg(__nu, __x);

	  const auto __factor
	    = __gnu_cxx::numbers::__root_pi_div_2_v<_Tp>
	      / std::sqrt(_Cmplx(__x));

	  const auto __j_n = __factor * _Bess.__J_value;
	  const auto __jp_n = __factor * _Bess.__J_deriv
			    - __j_n / (_Tp{2} * __x);
	  const auto __n_n = __factor * _Bess.__N_value;
	  const auto __np_n = __factor * _Bess.__N_deriv
			    - __n_n / (_Tp{2} * __x);

	  return __bess_t{__n, __x, __j_n, __jp_n, __n_n, __np_n};
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
      else if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__x == _Tp{0})
	{
	  if (__n == 0)
	    return _Tp{1};
	  else
	    return _Tp{0};
	}
      else
	return __sph_bessel_jn(__n, __x).__j_value;
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
      else if (std::isnan(__x))
	return __gnu_cxx::__quiet_NaN(__x);
      else if (__x == _Tp{0})
	return -__gnu_cxx::__infinity(__x);
      else
	return __sph_bessel_jn(__n, __x).__n_value;
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
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx _S_i{0, 1};
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__x);
      if (std::isnan(__x))
	return _Cmplx{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  const auto _Bess = __sph_bessel_jn_neg_arg(__n, __x);
	  return _Bess.__j_value + _S_i * _Bess.__n_value;
	}
      else
	{
	  const auto _Bess = __sph_bessel_jn(__n, __x);
	  return _Cmplx{_Bess.__j_value, _Bess.__n_value};
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
      using _Cmplx = std::complex<_Tp>;
      constexpr _Cmplx _S_i{0, 1};
      const auto _S_nan = __gnu_cxx::__quiet_NaN(__x);
      if (std::isnan(__x))
	return _Cmplx{_S_nan, _S_nan};
      else if (__x < _Tp{0})
	{
	  const auto _Bess = __sph_bessel_jn_neg_arg(__n, __x);
	  return _Bess.__j_value - _S_i * _Bess.__n_value;
	}
      else
	{
	  const auto _Bess = __sph_bessel_jn(__n, __x);
	  return _Cmplx{_Bess.__j_value, -_Bess.__n_value};
	}
    }

} // namespace __detail

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace std

#endif // _GLIBCXX_BITS_SF_BESSEL_TCC
