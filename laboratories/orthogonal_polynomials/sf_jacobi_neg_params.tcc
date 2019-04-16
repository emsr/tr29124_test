// Special functions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

#ifndef _LABS_SF_JACOBI_TCC
#define _LABS_SF_JACOBI_TCC 1

#include <ext/math_const.h>
#include <ext/math_util.h>

namespace lab
{

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
   * This works for @f$ \alpha,\beta > -1 @f$
   *
   * @tparam  _Tp  The real type of the argument
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first order parameter of the Jacobi polynomial
   * @param[in]  beta1  The second order parameter of the Jacobi polynomial
   * @param[in]  x  The argument
   */
  template<typename _Tp>
    __gnu_cxx::__jacobi_t<_Tp>
    __jacobi_recur(unsigned int __n, _Tp __alpha1, _Tp __beta1, _Tp __x)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__x);

      if (std::isnan(__alpha1) || std::isnan(__beta1) || std::isnan(__x))
	return {__n, __alpha1, __beta1, __x, _S_NaN, _S_NaN, _S_NaN};

      const auto __apb = __alpha1 + __beta1;
      auto __m = int(__n);
      if (const auto __pint = __gnu_cxx::__fp_is_integer(__n + 1 + __apb);
	  __pint && __pint() <= 0 && __m + 1 > -__pint())
	__m = -__pint();

      auto __P_nm2 = _Tp{1};
      if (__n == 0)
	return {__n, __alpha1, __beta1, __x, __P_nm2, _Tp{0}, _Tp{0}};

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
      for (auto __k = 3u; __k <= __n; ++__k )
	{
	  const auto __apbpk = __apb + _Tp(__k);
	  const auto __apbp2k = __apbpk + _Tp(__k);
	  const auto __apbp2km1 = __apbp2k - _Tp{1};
	  const auto __apbp2km2 = __apbp2km1 - _Tp{1};
	  const auto __d = _Tp{2} * __k * __apbpk * __apbp2km2;
	  const auto __a = __apbp2km2 * __apbp2km1 * __apbp2k;
	  const auto __b = __apbp2km1 * __a2mb2;
	  const auto __c = _Tp{2} * (__alpha1 + _Tp(__k - 1))
				  * (__beta1 + _Tp(__k - 1)) * __apbp2k;
	  if (__d == _Tp{0})
{
std::cerr << '\n' << '*' << ' ' << __n << ' ' << __alpha1 << ' ' << __beta1 << ' ' << __m << ' ' << __k << '\n';
break;
	//    std::__throw_runtime_error(__N("__jacobi_recur: "
	//				   "Failure in recursion"));
}
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
   * @f$ P_n^{(\alpha,\beta)}(x) @f$.
   * Thias works for @f$ \alpha, \beta > -1 @f$.
   *
   * @tparam  _Tp  The real type of the parameters
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first order parameter of the Jacobi polynomial
   * @param[in]  beta1  The second order parameter of the Jacobi polynomial
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
	      const auto __an = __alpha1 / __n;
	      const auto __bn = __beta1 / __n;
	      const auto __r1 = (1.0 + __alpha1) * (2.78 / (4.0 + __n * __n)
			+ 0.768 * __an / __n);
	      const auto __r2 = 1.0 + 1.48 * __an + 0.96 * __bn
			      + 0.452 * __an * __an + 0.83 * __an * __bn;
	      __z = 1.0 - __r1 / __r2;
	    }
	  else if (__i == 2)
	    {
	      const auto __r1 = (4.1 + __alpha1)
			/ ((1.0 + __alpha1) * (1.0 + 0.156 * __alpha1));
	      const auto __r2 = 1.0
			+ 0.06 * (__n - 8.0) * (1.0 + 0.12 * __alpha1) / __n;
	      const auto __r3 = 1.0
			+ 0.012 * __beta1 * (1.0 + 0.25 * std::abs(__alpha1))
			/ __n;
	      __z -= (1.0 - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == 3)
	    {
	      const auto __r1 = (1.67 + 0.28 * __alpha1)
			      / (1.0 + 0.37 * __alpha1);
	      const auto __r2 = 1.0 + 0.22 * (__n - 8.0) / __n;
	      const auto __r3 = 1.0 + 8.0 * __beta1
				/ ((6.28 + __beta1) * __n * __n);
	      __z -= (__pt[0].__point - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n - 1)
	    {
	      const auto __r1 = (1.0 + 0.235 * __beta1)
			      / (0.766 + 0.119 * __beta1);
	      const auto __r2 = 1.0 / (1.0 + 0.639 * (__n - 4.0)
						/ (1.0 + 0.71 * (__n - 4.0)));
	      const auto __r3 = 1.0 / (1.0 + 20.0 * __alpha1
				/ ((7.5 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 4].__point) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n)
	    {
	      const auto __r1 = (1.0 + 0.37 * __beta1)
			      / (1.67 + 0.28 * __beta1);
	      const auto __r2 = 1.0 / (1.0 + 0.22 * (__n - 8.0) / __n);
	      const auto __r3 = 1.0 / (1.0 + 8.0 * __alpha1
				/ ((6.28 + __alpha1) * __n * __n));
	      __z += (__z - __pt[__n - 3].__point) * __r1 * __r2 * __r3;
	    }
	  else
	    {
	      __z = 3.0 * __pt[__i - 2].__point
		  - 3.0 * __pt[__i - 3].__point + __pt[__i - 4].__point;
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
	  __pt[__i - 1].__point = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * @f$ n and m @f$, and real radial argument @f$ \rho @f$
   * is a polynomial of degree @f$ m + 2n @f$ in @f$ \rho @f$.
   *
   * The radial polynomials are defined by 
   * @f[
   *     R_n^m(\rho) = \sum_{k=0}^{\frac{n-m}{2}}
   *       \frac{(-1)^k(n-k)!}{k!(\frac{n+m}{2}-k)!(\frac{n-m}{2}-k)!}
   *       \rho^{n-2k}
   * @f]
   * for @f$ n - m @f$ even and identically 0 for @f$ n - m @f$ odd.
   *
   * The radial polynomials are related to the Jacobi polynomials:
   * @f[
   *    R_n^m(\rho) = (-1)^n x^m P_n^{(m,0)}(1-2\rho^2)
   * @f]
   * for @f$ 0 <= \rho <= 1 @f$
   *
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
   *      Cambridge University Press, 1999, pp 523-525 and 905-910.
   *
   * @see Zernike Polynomials: Evaluation, Quadrature, and Interpolation
   *      Philip Greengard, Kirill Serkh,
   *      Technical Report YALEU/DCS/TR-1539, February 20, 2018
   *
   * @tparam _Tp The real type of the radial coordinate
   * @param __n The non-negative degree.
   * @param __m The non-negative azimuthal order
   * @param __rho The radial argument
   */
  template<typename _Tp>
    _Tp
    __radial_jacobi(unsigned int __n, unsigned int __m, _Tp __rho)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__rho);

      if (std::isnan(__rho))
	return _S_NaN;

      const int __nmm = int(__n) - int(__m);
      if (__nmm < 0)
	return _Tp{0}; // FIXME: Is this true?
      else if (__nmm % 2 == 1)
	return _Tp{0};
      else
	{
	  auto __k = __nmm / 2;
	  return (__k % 2 == 0 ? +1 : -1) * std::pow(__rho, __m)
	       * __jacobi_recur(__k, _Tp(__m), _Tp{0},
				_Tp{1} - _Tp{2} * __rho * __rho).__P_n;
	}
    }

  /**
   * Return a vector containing the zeros of the radial Jacobi polynomial
   * @f$ P_n^{(\alpha,\beta)}(1 - 2\rho^2) @f$.
   *
   * @tparam _Tp The real type of the radial coordinate
   * @param[in]  n  The order of the Jacobi polynomial
   * @param[in]  alpha1  The first parameter of the Jacobi polynomial
   * @param[in]  beta1  The second parameter of the Jacobi polynomial
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __radial_jacobi_zeros(unsigned int __n, unsigned int __m)
    {
      const int __nmm = int(__n) - int(__m);
      if (__nmm % 2 != 0)
	return std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>();
      else
	{
	  auto __k = (int(__n) - int(__m)) / 2;
          auto __roots = __jacobi_zeros(__k, _Tp(__m), _Tp{0});
	  for (auto& __z : __roots)
	    {
	      __z.__point = std::sqrt((_Tp{1} - __z.__point) / _Tp{2});
	      // weight?
	    }
	  return __roots;
	}
    }

  /**
   * Return the Zernike polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative integral degree @f$ n @f$, signed integral order
   * @f$ m @f$, and real radial argument @f$ \rho @f$ and azimuthal angle
   * @f$ \phi @f$.
   *
   * The even Zernike polynomials are defined by:
   * @f[
   *    Z_n^m(\rho,\phi) = R_n^m(\rho)\cos(m\phi)
   * @f]
   * and the odd Zernike polynomials are defined by:
   * @f[
   *    Z_n^{-m}(\rho,\phi) = R_n^m(\rho)\sin(m\phi)
   * @f]
   * for non-negative degree @f$ m @f$ and @f$ m <= n @f$
   * and where @f$ R_n^m(\rho) @f$ is the radial polynomial
   * (@see __radial_jacobi).
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
    __gnu_cxx::fp_promote_t<_Tp>
    __zernike(unsigned int __n, int __m, _Tp __rho, _Tp __phi)
    {
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(__rho);

      if (std::isnan(__rho) || std::isnan(__phi))
	return _S_NaN;
      else
	return __radial_jacobi(__n, std::abs(__m), __rho)
	     * (__m >= 0 ? std::cos(__m * __phi) : std::sin(__m * __phi));
    }

} // namespace lab

#endif // _LABS_SF_JACOBI_TCC
