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

/** @file emsr/sf_jacobi.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef _LABS_SF_JACOBI_TCC
#define _LABS_SF_JACOBI_TCC 1

#include <emsr/math_constants.h>
#include <emsr/math_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/specfun_state.h>
#include <emsr/quadrature_point.h>
#include <emsr/special_functions.h>
#include <emsr/fp_type_util.h>

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
    emsr::jacobi_t<_Tp>
    jacobi_recur(unsigned int n, _Tp alpha1, _Tp beta1, _Tp x)
    {
      const auto _S_NaN = emsr::quiet_NaN(x);

      if (std::isnan(alpha1) || std::isnan(beta1) || std::isnan(x))
	return {n, alpha1, beta1, x, _S_NaN, _S_NaN, _S_NaN};

      const auto apb = alpha1 + beta1;
      auto m = int(n);
      if (const auto pint = emsr::fp_is_integer(n + 1 + apb);
	  pint && pint() <= 0 && m + 1 > -pint())
	m = -pint();

      auto P_nm2 = _Tp{1};
      if (n == 0)
	return {n, alpha1, beta1, x, P_nm2, _Tp{0}, _Tp{0}};

      const auto amb = alpha1 - beta1;
      auto P_nm1 = (amb + (apb + _Tp{2}) * x) / _Tp{2};
      if (n == 1)
	return {n, alpha1, beta1, x, P_nm1, P_nm2, _Tp{0}};

      const auto a2mb2 = amb * apb;
      const auto bah = ((apb + _Tp{2}) + _Tp{2});
      const auto poo = (bah - _Tp{1});
      auto P_n = (((poo * a2mb2)
		     + ((poo - _Tp{1}) * poo * bah) * x)
		    * P_nm1 - (_Tp{2} * (alpha1 + _Tp{1})
			    * (beta1 + _Tp{1}) * bah) * P_nm2)
		 / (_Tp{4} * (apb + _Tp{2}) * (poo - _Tp{1}));
      for (auto k = 3u; k <= n; ++k )
	{
	  const auto apbpk = apb + _Tp(k);
	  const auto apbp2k = apbpk + _Tp(k);
	  const auto apbp2km1 = apbp2k - _Tp{1};
	  const auto apbp2km2 = apbp2km1 - _Tp{1};
	  const auto d = _Tp{2} * k * apbpk * apbp2km2;
	  const auto a = apbp2km2 * apbp2km1 * apbp2k;
	  const auto b = apbp2km1 * a2mb2;
	  const auto c = _Tp{2} * (alpha1 + _Tp(k - 1))
				  * (beta1 + _Tp(k - 1)) * apbp2k;
	  if (d == _Tp{0})
{
std::cerr << '\n' << '*' << ' ' << n << ' ' << alpha1 << ' ' << beta1 << ' ' << m << ' ' << k << '\n';
break;
	//    throw std::runtime_error("jacobi_recur: Failure in recursion");
}
	  P_nm2 = P_nm1;
	  P_nm1 = P_n;
	  P_n = ((b + a * x) * P_nm1 - c * P_nm2) / d;
	}
      //auto Pp_n = (n * (alpha1 - beta1 - apbp2k * x) * P_nm1
	//	   + _Tp{2} * (n + alpha1) * (n + beta1) * P_nm2)
	//	/ (apbp2k * (_Tp{1} - x * x));
      return {n, alpha1, beta1, x, P_n, P_nm1, P_nm2};
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
    std::vector<emsr::QuadraturePoint<_Tp>>
    jacobi_zeros(unsigned int n, _Tp alpha1, _Tp beta1)
    {
      const auto _S_eps = emsr::epsilon(alpha1);
      const unsigned int _S_maxit = 1000u;

      std::vector<emsr::QuadraturePoint<_Tp>> pt(n);

      _Tp z;
      _Tp w;
      for (auto i = 1u; i <= n; ++i)
	{
	  if (i == 1)
	    {
	      const auto an = alpha1 / n;
	      const auto bn = beta1 / n;
	      const auto r1 = (1.0 + alpha1) * (2.78 / (4.0 + n * n)
			+ 0.768 * an / n);
	      const auto r2 = 1.0 + 1.48 * an + 0.96 * bn
			      + 0.452 * an * an + 0.83 * an * bn;
	      z = 1.0 - r1 / r2;
	    }
	  else if (i == 2)
	    {
	      const auto r1 = (4.1 + alpha1)
			/ ((1.0 + alpha1) * (1.0 + 0.156 * alpha1));
	      const auto r2 = 1.0
			+ 0.06 * (n - 8.0) * (1.0 + 0.12 * alpha1) / n;
	      const auto r3 = 1.0
			+ 0.012 * beta1 * (1.0 + 0.25 * std::abs(alpha1))
			/ n;
	      z -= (1.0 - z) * r1 * r2 * r3;
	    }
	  else if (i == 3)
	    {
	      const auto r1 = (1.67 + 0.28 * alpha1)
			      / (1.0 + 0.37 * alpha1);
	      const auto r2 = 1.0 + 0.22 * (n - 8.0) / n;
	      const auto r3 = 1.0 + 8.0 * beta1
				/ ((6.28 + beta1) * n * n);
	      z -= (pt[0].point - z) * r1 * r2 * r3;
	    }
	  else if (i == n - 1)
	    {
	      const auto r1 = (1.0 + 0.235 * beta1)
			      / (0.766 + 0.119 * beta1);
	      const auto r2 = 1.0 / (1.0 + 0.639 * (n - 4.0)
						/ (1.0 + 0.71 * (n - 4.0)));
	      const auto r3 = 1.0 / (1.0 + 20.0 * alpha1
				/ ((7.5 + alpha1) * n * n));
	      z += (z - pt[n - 4].point) * r1 * r2 * r3;
	    }
	  else if (i == n)
	    {
	      const auto r1 = (1.0 + 0.37 * beta1)
			      / (1.67 + 0.28 * beta1);
	      const auto r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
	      const auto r3 = 1.0 / (1.0 + 8.0 * alpha1
				/ ((6.28 + alpha1) * n * n));
	      z += (z - pt[n - 3].point) * r1 * r2 * r3;
	    }
	  else
	    {
	      z = 3.0 * pt[i - 2].point
		  - 3.0 * pt[i - 3].point + pt[i - 4].point;
	    }

	  auto alphabeta = alpha1 + beta1;
	  for (auto its = 1u; its <= _S_maxit; ++its)
	    {
	      auto temp = _Tp{2} + alphabeta;
	      auto P1 = (alpha1 - beta1 + temp * z) / _Tp{2};
	      auto P2 = _Tp{1};
	      for (auto j = 2u; j <= n; ++j)
		{
		  auto P3 = P2;
		  P2 = P1;
		  temp = _Tp{2} * j + alphabeta;
		  auto a = _Tp{2} * j * (j + alphabeta)
			   * (temp - _Tp{2});
		  auto b = (temp - _Tp{1})
			   * (alpha1 * alpha1 - beta1 * beta1
				+ temp * (temp - _Tp{2}) * z);
		  auto c = _Tp{2} * (j - 1 + alpha1)
			   * (j - 1 + beta1) * temp;
		  P1 = (b * P2 - c * P3) / a;
		}
	      auto Pp = (n * (alpha1 - beta1 - temp * z) * P1
			   + _Tp{2} * (n + alpha1) * (n + beta1) * P2)
			/ (temp * (_Tp{1} - z * z));
	      auto z1 = z;
	      z = z1 - P1 / Pp;
	      if (std::abs(z - z1) <= _S_eps)
		{
		  w = std::exp(std::lgamma(alpha1 + _Tp(n))
			       + std::lgamma(beta1 + _Tp(n))
			       - std::lgamma(_Tp(n + 1))
			       - std::lgamma(_Tp(n + 1) + alphabeta))
		      * temp * std::pow(_Tp{2}, alphabeta) / (Pp * P2);
		  break;
		}
	      if (its > _S_maxit)
		throw std::logic_error("jacobi_zeros: Too many iterations");
	    }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
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
   * @param n The non-negative degree.
   * @param m The non-negative azimuthal order
   * @param rho The radial argument
   */
  template<typename _Tp>
    _Tp
    radial_jacobi(unsigned int n, unsigned int m, _Tp rho)
    {
      const auto _S_NaN = emsr::quiet_NaN(rho);

      if (std::isnan(rho))
	return _S_NaN;

      const int nmm = int(n) - int(m);
      if (nmm < 0)
	return _Tp{0}; // FIXME: Is this true?
      else if (nmm % 2 == 1)
	return _Tp{0};
      else
	{
	  auto k = nmm / 2;
	  return (k % 2 == 0 ? +1 : -1) * std::pow(rho, m)
	       * jacobi_recur(k, _Tp(m), _Tp{0},
				_Tp{1} - _Tp{2} * rho * rho).P_n;
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
    std::vector<emsr::QuadraturePoint<_Tp>>
    radial_jacobi_zeros(unsigned int n, unsigned int m)
    {
      const int nmm = int(n) - int(m);
      if (nmm % 2 != 0)
	return std::vector<emsr::QuadraturePoint<_Tp>>();
      else
	{
	  auto k = (int(n) - int(m)) / 2;
          auto roots = jacobi_zeros(k, _Tp(m), _Tp{0});
	  for (auto& z : roots)
	    {
	      z.point = std::sqrt((_Tp{1} - z.point) / _Tp{2});
	      // weight?
	    }
	  return roots;
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
   * (@see radial_jacobi).
   *
   * @see Principals of Optics, 7th edition, Max Born and Emil Wolf,
   * Cambridge University Press, 1999, pp 523-525 and 905-910.
   *
   * @tparam _Tp The real type of the radial coordinate and azimuthal angle
   * @param n The non-negative integral degree.
   * @param m The integral azimuthal order
   * @param rho The radial coordinate
   * @param phi The azimuthal angle
   */
  template<typename _Tp>
    emsr::fp_promote_t<_Tp>
    zernike(unsigned int n, int m, _Tp rho, _Tp phi)
    {
      const auto _S_NaN = emsr::quiet_NaN(rho);

      if (std::isnan(rho) || std::isnan(phi))
	return _S_NaN;
      else
	return radial_jacobi(n, std::abs(m), rho)
	     * (m >= 0 ? std::cos(m * phi) : std::sin(m * phi));
    }

} // namespace lab

#endif // _LABS_SF_JACOBI_TCC
