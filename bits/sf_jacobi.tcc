// Special functions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
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

/** @file bits/sf_jacobi.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _GLIBCXX_BITS_SF_JACOBI_TCC
#define _GLIBCXX_BITS_SF_JACOBI_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Compute the Jacobi polynomial by recursion on @c n:
   * @f[
   *   2 n(\alpha + \beta + n) (\alpha + \beta + 2n - 2)
   *         P^{(\alpha, \beta)}_{n}(x)
   *     = (\alpha + \beta + 2n - 1)
   *       ((\alpha^2 - \beta^2)
   *        + x(\alpha + \beta + 2n - 2)(\alpha + \beta + 2n))
   *         P^{(\alpha, \beta)}_{n-1}(x)
   *     - 2 (\alpha + n - 1)(\beta + n - 1)(\alpha + \beta + 2n)
   *         P^{(\alpha, \beta)}_{n-2}(x)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __poly_jacobi(unsigned int __n, _Tp __alpha, _Tp __beta, _Tp __x)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();

      if (__isnan(__alpha) || __isnan(__beta) || __isnan(__x))
	return _S_NaN;

      auto _Pm2 = _Tp{1};
      if (__n == 0)
	return _Pm2;

      auto __apb = __alpha + __beta;
      auto __amb = __alpha - __beta;
      auto _Pm1 = (__amb + (_Tp{2} + __apb) * __x) / _Tp{2};
      if (__n == 1)
	return _Pm1;

      auto _Pm0 = _Tp{0};
      auto __a2mb2 = __amb * __apb;
      for (auto __k = 2; __k <= __n; ++__k )
	{
	  auto __apbpk = __apb + _Tp(__k);
	  auto __apbp2k = __apbpk + _Tp(__k);
	  auto __apbp2km1 = __apbp2k - _Tp{1};
	  auto __apbp2km2 = __apbp2km1 - _Tp{1};
	  auto __d = _Tp{2} * __k * __apbpk * __apbp2km2;
	  auto __a = __apbp2km2 * __apbp2km1 * __apbp2k;
	  auto __b = __apbp2km1 * __a2mb2;
	  auto __c = _Tp{2} * (__alpha + _Tp(__k - 1))
			    * (__beta + _Tp(__k - 1)) * __apbp2k;
	  if (__d == _Tp{0})
	    std::__throw_runtime_error(__N("__poly_jacobi: "
					   "Failure in recursion"));
	  _Pm0 = ((__b + __a * __x) * _Pm1 - __c * _Pm2) / __d;
	  _Pm2 = _Pm1;
	  _Pm1 = _Pm0;
	}
      return _Pm0;
    }
/*

  P^{(a,b)}_{k}(z) =

    (a+b+2k-1)((a-b)(a+b) + z(a+b+2k-2)(a+b+2k))
    -------------------------------------------- P^{(a,b)}_{k-1}(z)
             2 k(a+b+k)(a+b+2k-2)

       2(a+k-1)(b+k-1)(a+b+2k)
    -  ----------------------- P^{(a,b)}_{k-2}(z)
	 2 k(a+b+k)(a+b+2k-2)

*/

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @f$ n @f$, order @f$ m <= n @f$, and real radial
   * argument @f$ \rho @f$.
   *
   * The radial polynomials are defined by 
   * @f[
   *     R_n^m(\rho) = \sum_{k=0}^{\frac{n-m}{2}}
   *       \frac{(-1)^k(n-k)!}{k!(\frac{n+m}{2}-k)!(\frac{n-m}{2}-k)!}
   *       \rho^{n-2k}
   * @f]
   * for @f$ n - m @f$ even and identically 0 for @f$ n - m @f$ odd.
   * The radial polynomials can be related to the Zernike polynomials:
   * @f[
   *    Z_n^m(\rho,\phi) = R_n^m(\rho) \cos(m\phi)
   * @f]
   * @f[
   *    Z_n^{-m}(\rho,\phi) = R_n^m(\rho) \sin(m\phi)
   * @f]
   * for non-negative  @f$ m, n @f$.
   * @see zernike for details on the Zernike polynomials.
   *
   * @see Principals of Optics, 7th edition, Max Born and Emil Wolf,
   * Cambridge University Press, 1999, pp 523-525 and 905-910.
   *
   * @tparam _Tp The real type of the radial coordinate
   * @param __n The non-negative degree.
   * @param __m The non-negative azimuthal order
   * @param __rho The radial argument
   */
  template<typename _Tp>
    _Tp
    __poly_radial_jacobi(unsigned int __n, unsigned int __m, _Tp __rho)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();

      if (__isnan(__rho))
	return _S_NaN;

      if (__m > __n)
	std::__throw_range_error(__N("poly_radial_jacobi: order > degree"));
      else if ((__n - __m) % 2 == 1)
	return _Tp{0};
      else
	{
	  auto __k = (__n - __m) / 2;
	  return (__k % 2 == 0 ? +1 : -1) * std::pow(__rho, __m)
	       * __poly_jacobi(__k, _Tp(__m), _Tp{0},
			       _Tp{1} - _Tp{2} * __rho * __rho);
	}
    }

  /**
   * Return the Zernicke polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative integral degree @f$ n @f$, signed integral order
   * @f$ m @f$, and real radial argument @f$ \rho @f$ and azimuthal angle
   * @f$ \phi @f$.
   *
   * The even Zernicke polynomials are defined by:
   * @f[
   *    Z_n^m(\rho,\phi) = R_n^m(\rho)\cos(m\phi)
   * @f]
   * and the odd Zernicke polynomials are defined by:
   * @f[
   *    Z_n^{-m}(\rho,\phi) = R_n^m(\rho)\sin(m\phi)
   * @f]
   * for non-negative degree @f$ m @f$ and @f$ m <= n @f$
   * and where @f$ R_n^m(\rho) @f$ is the radial polynomial
   * (@see __poly_radial_jacobi).
   *
   * @see Principals of Optics, 7th edition, Max Born and Emil Wolf,
   * Cambridge University Press, 1999, pp 523-525 and 905-910.
   *
   * @tparam _Tp The real type of the radial coordinate and azimuthal angle
   * @param __n The non-negative integral degree.
   * @param __m The integral azimuthal order
   * @param __rho The radial coordinate
   * @param __phi The azimuthal angle
   */
  template<typename _Tp>
    __gnu_cxx::__promote_fp_t<_Tp>
    __zernike(unsigned int __n, int __m, _Tp __rho, _Tp __phi)
    {
      constexpr auto _S_NaN = __gnu_cxx::__quiet_NaN<_Tp>();

      if (__isnan(__rho) || __isnan(__phi))
	return _S_NaN;
      else
        return __poly_radial_jacobi(__n, std::abs(__m), __rho)
	     * (__m >= 0 ? std::cos(__m * __phi) : std::sin(__m * __phi));
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_JACOBI_TCC
