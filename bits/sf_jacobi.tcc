// Special functions -*- C++ -*-

// Copyright (C) 2016-2017 Free Software Foundation, Inc.
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
    __gnu_cxx::__jacobi_t<_Tp>
    __jacobi_recur(unsigned int __n, _Tp __alpha1, _Tp __beta1, _Tp __x)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);

      if (__isnan(__alpha1) || __isnan(__beta1) || __isnan(__x))
	return {__n, __alpha1, __beta1, _S_NaN, _S_NaN, _S_NaN, _S_NaN};

      auto __P_nm2 = _Tp{1};
      if (__n == 0)
	return {__n, __alpha1, __beta1, __x, __P_nm2, _Tp{0}, _Tp{0}};

      const auto __apb = __alpha1 + __beta1;
      const auto __amb = __alpha1 - __beta1;
      auto __P_nm1 = (__amb + (__apb + _Tp{2}) * __x) / _Tp{2};
      if (__n == 1)
	return {__n, __alpha1, __beta1, __x, __P_nm1, __P_nm2, _Tp{0}};

      const auto __a2mb2 = __amb * __apb;
      const auto __bah = ((__apb + _Tp{2}) + _Tp{2});
      const auto __poo = (__bah - _Tp{1});
      auto __P_n = (((__poo * __a2mb2)
		     + ((__poo - _Tp{1}) * __poo * __bah) * __x)
		    * __P_nm1 - (_Tp{2} * (__alpha1 + _Tp{1})
			    * (__beta1 + _Tp{1}) * __bah) * __P_nm2)
		 / (_Tp{4} * (__apb + _Tp{2}) * (__poo - _Tp{1}));
      for (auto __k = 3; __k <= __n; ++__k )
	{
	  auto __apbpk = __apb + _Tp(__k);
	  auto __apbp2k = __apbpk + _Tp(__k);
	  auto __apbp2km1 = __apbp2k - _Tp{1};
	  auto __apbp2km2 = __apbp2km1 - _Tp{1};
	  auto __d = _Tp{2} * __k * __apbpk * __apbp2km2;
	  auto __a = __apbp2km2 * __apbp2km1 * __apbp2k;
	  auto __b = __apbp2km1 * __a2mb2;
	  auto __c = _Tp{2} * (__alpha1 + _Tp(__k - 1))
			    * (__beta1 + _Tp(__k - 1)) * __apbp2k;
	  if (__d == _Tp{0})
	    std::__throw_runtime_error(__N("__jacobi_recur: "
					   "Failure in recursion"));
	  __P_nm2 = __P_nm1;
	  __P_nm1 = __P_n;
	  __P_n = ((__b + __a * __x) * __P_nm1 - __c * __P_nm2) / __d;
	}
      //auto __Pp_n = (__n * (__alpha1 - __beta1 - __apbp2k * __x) * __P_nm1
	//	   + _Tp{2} * (__n + __alpha1) * (__n + __beta1) * __P_nm2)
	//	/ (__apbp2k * (_Tp{1} - __x * __x));
      return {__n, __alpha1, __beta1, __x, __P_n, __P_nm1, __P_nm2};
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
   * Return a vector containing the zeros of the Jacobi polynomial
   * @f$ P_n^{(\alpha,\beta)}@f$.
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __jacobi_zeros(unsigned int __n, _Tp __alpha1, _Tp __beta1)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha1);
      const unsigned int _S_maxit = 1000u;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);

      _Tp __z;
      _Tp __w;
      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  if (__i == 1)
	    {
	      auto __an = __alpha1 / __n;
	      auto __bn = __beta1 / __n;
	      auto __r1 = (1.0 + __alpha1) * (2.78 / (4.0 + __n * __n)
			+ 0.768 * __an / __n);
	      auto __r2 = 1.0 + 1.48 * __an + 0.96 * __bn
			+ 0.452 * __an * __an + 0.83 * __an * __bn;
	      __z = 1.0 - __r1 / __r2;
	    }
	  else if (__i == 2)
	    {
	      auto __r1 = (4.1 + __alpha1)
			/ ((1.0 + __alpha1) * (1.0 + 0.156 * __alpha1));
	      auto __r2 = 1.0
			+ 0.06 * (__n - 8.0) * (1.0 + 0.12 * __alpha1) / __n;
	      auto __r3 = 1.0
		    + 0.012 * __beta1 * (1.0 + 0.25 * std::abs(__alpha1)) / __n;
	      __z -= (1.0 - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == 3)
	    {
	      auto __r1 = (1.67 + 0.28 * __alpha1) / (1.0 + 0.37 * __alpha1);
	      auto __r2 = 1.0 + 0.22 * (__n - 8.0) / __n;
	      auto __r3 = 1.0 + 8.0 * __beta1 / ((6.28 + __beta1) * __n * __n);
	      __z -= (__pt[0].__zero - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n - 1)
	    {
	      auto __r1 = (1.0 + 0.235 * __beta1) / (0.766 + 0.119 * __beta1);
	      auto __r2 = 1.0 / (1.0 + 0.639 * (__n - 4.0)
						/ (1.0 + 0.71 * (__n - 4.0)));
	      auto __r3 = 1.0 / (1.0 + 20.0 * __alpha1
				/ ((7.5 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 4].__zero) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n)
	    {
	      auto __r1 = (1.0 + 0.37 * __beta1) / (1.67 + 0.28 * __beta1);
	      auto __r2 = 1.0 / (1.0 + 0.22 * (__n - 8.0) / __n);
	      auto __r3 = 1.0 / (1.0 + 8.0 * __alpha1
				 / ((6.28 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 3].__zero) * __r1 * __r2 * __r3;
	    }
	  else
	    {
	      __z = 3.0 * __pt[__i - 2].__zero
		  - 3.0 * __pt[__i - 3].__zero + __pt[__i - 4].__zero;
	    }

	  auto __alphabeta = __alpha1 + __beta1;
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __temp = _Tp{2} + __alphabeta;
	      auto __P1 = (__alpha1 - __beta1 + __temp * __z) / _Tp{2};
	      auto __P2 = _Tp{1};
	      for (auto __j = 2u; __j <= __n; ++__j)
		{
		  auto __P3 = __P2;
		  __P2 = __P1;
		  __temp = _Tp{2} * __j + __alphabeta;
		  auto __a = _Tp{2} * __j * (__j + __alphabeta)
			   * (__temp - _Tp{2});
		  auto __b = (__temp - _Tp{1})
			   * (__alpha1 * __alpha1 - __beta1 * __beta1
				+ __temp * (__temp - _Tp{2}) * __z);
		  auto __c = _Tp{2} * (__j - 1 + __alpha1)
			   * (__j - 1 + __beta1) * __temp;
		  __P1 = (__b * __P2 - __c * __P3) / __a;
		}
	      auto __Pp = (__n * (__alpha1 - __beta1 - __temp * __z) * __P1
			   + _Tp{2} * (__n + __alpha1) * (__n + __beta1) * __P2)
			/ (__temp * (_Tp{1} - __z * __z));
	      auto __z1 = __z;
	      __z = __z1 - __P1 / __Pp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = std::exp(std::lgamma(__alpha1 + _Tp(__n))
			       + std::lgamma(__beta1 + _Tp(__n))
			       - std::lgamma(_Tp(__n + 1))
			       - std::lgamma(_Tp(__n + 1) + __alphabeta))
		      * __temp * std::pow(_Tp{2}, __alphabeta) / (__Pp * __P2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__jacobi_zeros: Too many iterations");
	    }
	  __pt[__i - 1].__zero = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

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
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__rho);

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
	       * __jacobi_recur(__k, _Tp(__m), _Tp{0},
				_Tp{1} - _Tp{2} * __rho * __rho).__P_n;
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
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__rho);

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
