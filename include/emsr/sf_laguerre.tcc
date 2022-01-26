
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

/** @file bits/sf_laguerre.tcc
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
//     Section 13, pp. 509-510, Section 22 pp. 773-802
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl

#ifndef SF_LAGUERRE_TCC
#define SF_LAGUERRE_TCC 1

#include <stdexcept>
#include <vector>

#include <emsr/specfun_state.h>
#include <emsr/quadrature_point.h>
#include <emsr/sf_trig.h> // sin_pi, cos_pi, polar_pi

namespace emsr
{
namespace detail
{

  /**
   * @brief This routine returns the associated Laguerre polynomial
   *        of degree @f$ n @f$, order @f$ \alpha > -1 @f$ for large n.
   * Abramowitz & Stegun, 13.5.21
   *
   * @tparam _Tpa The type of the order.
   * @tparam _Tp The type of the parameter.
   * @param n The degree of the Laguerre function.
   * @param alpha1 The order of the Laguerre function.
   * @param x The argument of the Laguerre function.
   * @return The value of the Laguerre function of degree n,
   *         order @f$ \alpha @f$, and argument x.
   *
   *  This is from the GNU Scientific Library.
   */
  template<typename _Tpa, typename _Tp>
    _Tp
    laguerre_large_n(unsigned n, _Tpa alpha1, _Tp x)
    {
      const auto a = -_Tp(n);
      const auto b = _Tp(alpha1) + _Tp{1};
      const auto eta = _Tp{2} * b - _Tp{4} * a;
      const auto cos2th = x / eta;
      const auto sin2th = _Tp{1} - cos2th;
      const auto th = std::acos(std::sqrt(cos2th));
      const auto pre_h = (emsr::pi_v<_Tp> / _Tp{2})
			* (emsr::pi_v<_Tp> / _Tp{2})
			* eta * eta * cos2th * sin2th;

      const auto lg_b = log_gamma(_Tp(n) + b);
      const auto lnfact = log_gamma(_Tp(n + 1));

      const auto pre_term1 = _Tp{0.5L} * (_Tp{1} - b)
			     * std::log(_Tp{0.25L} * x * eta);
      const auto pre_term2 = _Tp{0.25L} * std::log(pre_h);
      const auto lnpre = lg_b - lnfact + _Tp{0.5L} * x
			 + pre_term1 - pre_term2;
      const auto ser_term1 = emsr::sin_pi(a);
      const auto ser_term2 = std::sin(_Tp{0.25L} * eta
			     * (_Tp{2} * th - std::sin(_Tp{2} * th))
			      + emsr::pi_v<_Tp> / _Tp{4});
      const auto ser = ser_term1 + ser_term2;

      return std::exp(lnpre) * ser;
    }


  /**
   * @brief Evaluate the polynomial based on the confluent hypergeometric
   *        function in a safe way, with no restriction on the arguments.
   *
   * The associated Laguerre function is defined by
   * @f[
   *    L_n^{(\alpha)}(x) = \frac{(\alpha + 1)_n}{n!}
   * 		       {}_1F_1(-n; \alpha + 1; x)
   * @f]
   * where @f$ (\alpha)_n @f$ is the Pochhammer symbol and
   * @f$ {}_1F_1(a; c; x) @f$ is the confluent hypergeometric function.
   *
   * This function assumes x != 0.
   *
   * This is from the GNU Scientific Library.
   *
   * @tparam _Tpa The type of the order.
   * @tparam _Tp The type of the parameter.
   * @param n The degree of the Laguerre function.
   * @param alpha1 The order of the Laguerre function.
   * @param x The argument of the Laguerre function.
   * @return The value of the Laguerre function of degree n,
   * 	     order @f$ \alpha @f$, and argument x.
   */
  template<typename _Tpa, typename _Tp>
    _Tp
    laguerre_hyperg(unsigned int n, _Tpa alpha1, _Tp x)
    {
      const auto b = _Tp(alpha1) + _Tp{1};
      const auto mx = -x;
      const auto tc_sgn = (x < _Tp{0} ? _Tp{1}
			 : ((n % 2 == 1) ? -_Tp{1} : _Tp{1}));
      // Get |x|^n/n!
      auto tc = _Tp{1};
      const auto ax = std::abs(x);
      for (unsigned int k = 1; k <= n; ++k)
	tc *= (ax / k);

      auto term = tc * tc_sgn;
      auto sum = term;
      for (int k = int(n) - 1; k >= 0; --k)
	{
	  term *= ((b + _Tp(k)) / _Tp(int(n) - k))
		  * _Tp(k + 1) / mx;
	  sum += term;
	}

      return sum;
    }


  /**
   * @brief This routine returns the associated Laguerre polynomial
   * 	    of degree @c n, order @c @f$ \alpha @f$: @f$ L_n^{(\alpha)}(x) @f$
   * 	    by recursion.
   *
   * The associated Laguerre function is defined by
   * @f[
   *   L_n^{(\alpha)}(x) = \frac{(\alpha + 1)_n}{n!}
   *                  {}_1F_1(-n; \alpha + 1; x)
   * @f]
   * where @f$ (\alpha)_n @f$ is the Pochhammer symbol and
   * @f$ {}_1F_1(a; c; x) @f$ is the confluent hypergeometric function.
   *
   * The associated Laguerre polynomial is defined for integral
   * order @f$ \alpha = m @f$ by:
   * @f[
   *    L_n^{(m)}(x) = (-1)^m \frac{d^m}{dx^m} L_{n + m}(x)
   * @f]
   * where the Laguerre polynomial is defined by:
   * @f[
   *    L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @tparam _Tpa The type of the order.
   * @tparam _Tp The type of the parameter.
   * @param n The degree of the Laguerre function.
   * @param alpha1 The order of the Laguerre function.
   * @param x The argument of the Laguerre function.
   * @return The value of the Laguerre function of order n,
   * 	     degree @f$ \alpha @f$, and argument x.
   */
  template<typename _Tpa, typename _Tp>
    emsr::laguerre_t<_Tpa, _Tp>
    laguerre_recur(unsigned int n, _Tpa alpha1, _Tp x)
    {
      // Compute L_0.
      auto L_0 = _Tp{1};
      if  (n == 0)
	return {n, alpha1, x, L_0, _Tp{0}, _Tp{0}};

      // Compute L_1^{(alpha)}.
      auto L_1 = -x + _Tp{1} + _Tp(alpha1);
      if  (n == 1)
	return {n, alpha1, x, L_1, L_0, _Tp{0}};

      // Compute L_n^{(alpha)} by recursion on n.
      auto L_nm2 = L_0;
      auto L_nm1 = L_1;
      auto L_n = (_Tp{3} + _Tp(alpha1) - x) * L_nm1 / _Tp{2}
		  - (_Tp{1} + _Tp(alpha1)) * L_nm2 / _Tp{2};
      for  (unsigned int nn = 3; nn <= n; ++nn)
	{
	    L_nm2 = L_nm1;
	    L_nm1 = L_n;
	    L_n = (_Tp(2 * nn - 1) + _Tp(alpha1) - x)
		  * L_nm1 / _Tp(nn)
		  - (_Tp(nn - 1) + _Tp(alpha1)) * L_nm2 / _Tp(nn);
	}

      // Derivative.
      //auto Lp_n = (_Tp(n) * L_nm1 - _Tp(n + alpha1) * L_nm2) / x;
      return {n, alpha1, x, L_n, L_nm1, L_nm2};
    }

  /**
   * Return an array of abscissae and weights for the Gauss-Laguerre rule.
   */
  template<typename _Tp>
    std::vector<emsr::QuadraturePoint<_Tp>>
    laguerre_zeros(unsigned int n, _Tp alpha1)
    {
      const auto s_eps = emsr::epsilon(alpha1);
      const unsigned int s_maxit = 1000;

      std::vector<emsr::QuadraturePoint<_Tp>> pt(n);

      for (auto i = 1u; i <= n; ++i)
	{
	  auto z = _Tp{0};
	  auto w = _Tp{0};
	  // Clever approximations for roots.
	  if (i == 1)
	    z += (1.0 + alpha1)
		 * (3.0 + 0.92 * alpha1) / (1.0 + 2.4 * n + 1.8 * alpha1);
	  else if (i == 2)
	    z += (15.0 + 6.25 * alpha1) / (1.0 + 2.5 * n + 0.9 * alpha1);
	  else
	    {
	      auto ai = i - 2;
	      z += ((1.0 + 2.55 * ai) / (1.9 * ai)
		     + 1.26 * ai * alpha1 / (1.0 + 3.5 * ai))
		   * (z - pt[i - 3].point) / (1.0 + 0.3 * alpha1);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto L2 = _Tp{0};
	      auto L1 = _Tp{1};
	      for (auto j = 1u; j <= n; ++j)
		{
		  auto L3 = L2;
		  L2 = L1;
		  L1 = ((_Tp(2 * j - 1 + alpha1) - z) * L2
			- (_Tp(j - 1 + alpha1)) * L3) / _Tp(j);
		}
	      // Derivative.
	      auto Lp = (_Tp(n) * L1 - _Tp(n + alpha1) * L2) / z;
	      // Newton's rule for root.
	      auto z1 = z;
	      z = z1 - L1 / Lp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  auto exparg = std::lgamma(_Tp(alpha1 + n))
				- std::lgamma(_Tp(n));
		  w = -std::exp(exparg) / (Lp * n * L2);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("laguerre_zeros: Too many iterations");
	   }
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}
    return pt;
  }


  /**
   * @brief This routine returns the associated Laguerre polynomial
   * 	    of degree n, order @f$ \alpha @f$: @f$ L_n^{(\alpha)}(x) @f$.
   *
   * The associated Laguerre function is defined by
   * @f[
   *    L_n^{(\alpha)}(x) = \frac{(\alpha + 1)_n}{n!}
   * 			{}_1F_1(-n; \alpha + 1; x)
   * @f]
   * where @f$ (\alpha)_n @f$ is the Pochhammer symbol and
   * @f$ {}_1F_1(a; c; x) @f$ is the confluent hypergeometric function.
   *
   * The associated Laguerre polynomial is defined for integral
   * order @f$ \alpha = m @f$ by:
   * @f[
   * 	 L_n^{(m)}(x) = (-1)^m \frac{d^m}{dx^m} L_{n + m}(x)
   * @f]
   * where the Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @tparam _Tpa The type of the order.
   * @tparam _Tp The type of the parameter.
   * @param n The degree of the Laguerre function.
   * @param alpha1 The order of the Laguerre function.
   * @param x The argument of the Laguerre function.
   * @return The value of the Laguerre function of degree n,
   *         order @f$ \alpha @f$, and argument x.
   */
  template<typename _Tpa, typename _Tp>
    _Tp
    laguerre(unsigned int n, _Tpa alpha1, _Tp x)
    {
      const unsigned int n_huge = 10000000;
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (n == 0)
	return _Tp{1};
      else if (n == 1)
	return _Tp{1} + _Tp(alpha1) - x;
      else if (x == _Tp{0})
	{
	  auto prod = _Tp(alpha1) + _Tp{1};
	  for (unsigned int k = 2; k <= n; ++k)
	    prod *= (_Tp(alpha1) + _Tp(k)) / _Tp(k);
	  return prod;
	}
      else if (n > n_huge && _Tp(alpha1) > -_Tp{1}
	    && x < _Tp{2} * (_Tp(alpha1) + _Tp{1}) + _Tp(4 * n))
	return laguerre_large_n(n, alpha1, x);
      else if (_Tp(alpha1) >= _Tp{0}
	   || (x > _Tp{0} && _Tp(alpha1) < -_Tp(n + 1)))
	return laguerre_recur(n, alpha1, x).L_n;
      else
	return laguerre_hyperg(n, alpha1, x);
    }


  /**
   * @brief This routine returns the associated Laguerre polynomial
   * 	    of degree n, order m: @f$ L_n^{(m)}(x) @f$.
   *
   * The associated Laguerre polynomial is defined for integral
   * order @f$ \alpha = m @f$ by:
   * @f[
   * 	 L_n^{(m)}(x) = (-1)^m \frac{d^m}{dx^m} L_{n + m}(x)
   * @f]
   * where the Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @tparam _Tpa The type of the order.
   * @tparam _Tp The type of the parameter
   * @param n The degree
   * @param alpha The order
   * @param x The argument
   * @return The value of the associated Laguerre polynomial of order n,
   *         degree m, and argument x.
   */
  template<typename _Tpa, typename _Tp>
    _Tp
    assoc_laguerre(unsigned int n, _Tpa alpha, _Tp x)
    { return laguerre<_Tpa, _Tp>(n, alpha, x); }


  /**
   * @brief This routine returns the Laguerre polynomial
   * 	    of degree n: @f$ L_n(x) @f$.
   *
   * The Laguerre polynomial is defined by:
   * @f[
   *    L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @param n The degree of the Laguerre polynomial.
   * @param x The argument of the Laguerre polynomial.
   * @return The value of the Laguerre polynomial of order n
   * 	     and argument x.
   */
  template<typename _Tp>
    _Tp
    laguerre(unsigned int n, _Tp x)
    { return laguerre<unsigned int, _Tp>(n, 0, x); }

} // namespace detail
} // namespace emsr

#endif // SF_LAGUERRE_TCC
