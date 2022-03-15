
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

/** @file bits/sf_hermite.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// Reference:
// (1) Handbook of Mathematical Functions,
//     Ed. Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, Section 22 pp. 773-802

#ifndef SF_HERMITE_TCC
#define SF_HERMITE_TCC 1

#include <stdexcept>

#include <emsr/sf_mod_bessel.h>

namespace emsr
{
namespace detail
{

  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of order n: @f$ H_n(x) @f$ by recursion on n.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * The Hermite polynomial has first and second derivatives:
   * @f[
   *    H'_n(x) = 2n H_{n-1}(x)
   * @f]
   * and
   * @f[
   *    H''_n(x) = 4n(n - 1) H_{n-2}(x)
   * @f]
   *
   * The Physicists Hermite polynomials have highest-order coefficient
   * @f$ 2^n @f$ and are orthogonal with respect to the weight function
   * @f[
   *   w(x) = e^{x^2}
   * @f]
   *
   * @param n The order of the Hermite polynomial.
   * @param x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename Tp>
    emsr::hermite_t<Tp>
    hermite_recur(unsigned int n, Tp x)
    {
      // Compute H_0.
      auto H_nm2 = Tp{1};
      if (n == 0)
	return {n, x, H_nm2, Tp{0}, Tp{0}};

      // Compute H_1.
      auto H_nm1 = Tp{2} * x;
      if (n == 1)
	return {n, x, H_nm1, H_nm2, Tp{0}};

      // Try to detect blowup.
      /// @todo Find the sign of Hermite blowup values.
      {
	constexpr auto s_inf = std::numeric_limits<Tp>::infinity();
	constexpr auto s_max = std::numeric_limits<Tp>::max();
	// Subtract little safety below.
	const auto arg = std::log2(s_max) / n - Tp{0.125};
	const auto xm = std::pow(Tp{2}, arg);
	const auto hinf = Tp(n & 1 ? -1 : +1) * s_inf;
	if (std::abs(x) > xm)
	  return {n, x, hinf, -hinf, hinf};
      }

      // Compute H_n.
      auto H_n = Tp{2} * (x * H_nm1 - H_nm2);
      for (unsigned int i = 3; i <= n; ++i)
	{
	  H_nm2 = H_nm1;
	  H_nm1 = H_n;
	  H_n = Tp{2} * (x * H_nm1 - Tp(i - 1) * H_nm2);
	}

      return {n, x, H_n, H_nm1, H_nm2};
    }

  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of large order n: @f$ H_n(x) @f$.  We assume here
   * 	    that x >= 0.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * @see "Asymptotic analysis of the Hermite polynomials
   * 	  from their differential-difference equation", 
   * 	  Diego Dominici, arXiv:math/0601078v1 [math.CA] 4 Jan 2006
   * @param n The order of the Hermite polynomial.
   * @param x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename Tp>
    Tp
    hermite_asymp(unsigned int n, Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_sqrt_2pi = emsr::sqrttau_v<Tp>;
      // x >= 0 in this routine.
      const auto xturn = std::sqrt(Tp(2 * n));
      if (std::abs(x - xturn) < Tp{0.05L} * xturn)
	{
	  // Transition region x ~ sqrt(2n).
	  const auto n_2 = Tp(n) / Tp{2};
	  const auto n6th = std::pow(Tp(n), Tp{1} / Tp{6});
	  const auto exparg = n * std::log(xturn) - Tp{3} * n_2
			      + xturn * x;
	  const auto airyarg = s_sqrt_2 * (x - xturn) * n6th;
	  auto Ai = emsr::airy_ai(airyarg);
	  return s_sqrt_2pi * n6th * std::exp(exparg) * Ai;
	}
      else if (x < xturn)
	{
	  // Oscillatory region |x| < sqrt(2n).
	  const auto theta = std::asin(x / xturn);
	  const auto __2theta = Tp{2} * theta;
	  const auto n_2 = Tp(n) / Tp{2};
	  const auto exparg = n_2 * (Tp{2} * std::log(xturn)
					- std::cos(__2theta));
	  const auto arg = theta / Tp{2}
		+ n_2 * (std::sin(__2theta) + __2theta - s_pi);
	  return std::sqrt(Tp{2} / std::cos(theta))
	       * std::exp(exparg) * std::cos(arg);
	}
      else
	{
	  // Exponential region |x| > sqrt(2n).
	  const auto sigma = std::sqrt((x - xturn) * (x + xturn));
	  const auto exparg = Tp{0.5L} * (x * (x - sigma) - n)
			     + n * std::log(sigma + x);
	  return std::exp(exparg)
	       * std::sqrt(Tp{0.5L} * (Tp{1} + x / sigma));
	}
    }


  /**
   * @brief This routine returns the Hermite polynomial
   * 	    of order n: @f$ H_n(x) @f$.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   * An explicit series formula is:
   * @f[
   *   H_n(x) = \sum_{k=0}^{m} \frac{(-1)^k}{k!(n-2k)!}(2x)^{n-2k}
   *     \mbox{ where } m = \left\lfloor{\frac{n}{2}}\right\rfloor
   * @f]
   *
   * The Hermite polynomial obeys a reflection formula:
   * @f[
   *   H_n(-x) = (-1)^n H_n(x)
   * @f]
   *
   * @param n The order of the Hermite polynomial.
   * @param x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename Tp>
    Tp
    hermite(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return (n % 2 == 1 ? -1 : +1) * hermite(n, -x);
      else if (n > 10000)
	return hermite_asymp(n, x);
      else
	return hermite_recur(n, x).H_n;
    }

  /**
   * @brief This routine returns the Probabilists Hermite polynomial
   * 	    of order n: @f$ He_n(x) @f$ by recursion on n.
   *
   * The Probabilists Hermite polynomial is defined by:
   * @f[
   *   He_n(x) = (-1)^n e^{x^2/2} \frac{d^n}{dx^n} e^{-x^2/2}
   * @f]
   * or
   * @f[
   *   He_n(x) = \frac{1}{2^{-n/2}}H_n\left(\frac{x}{\sqrt{2}}\right)
   * @f]
   * where @f$ H_n(x) @f$ is the Physicists Hermite function.
   *
   * The Probabilists Hermite polynomial has first and second derivatives:
   * @f[
   *    He'_n(x) = n He_{n-1}(x)
   * @f]
   * and
   * @f[
   *    He''_n(x) = n(n - 1) He_{n-2}(x)
   * @f]
   *
   * The Probabilists Hermite polynomial are monic and are orthogonal
   * with respect to the weight function
   * @f[
   *   w(x) = e^{x^2/2}
   * @f]
   *
   * @param n The order of the Hermite polynomial.
   * @param x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename Tp>
    emsr::hermite_he_t<Tp>
    prob_hermite_recur(unsigned int n, Tp x)
    {
      // Compute He_0.
      auto He_nm2 = Tp{1};
      if (n == 0)
	return {n, x, He_nm2, Tp{0}, Tp{0}};

      // Compute He_1.
      auto He_nm1 = x;
      if (n == 1)
	return {n, x, He_nm1, He_nm2, Tp{0}};

      // Compute He_n.
      auto He_n = x * He_nm1 - He_nm2;
      for (unsigned int i = 3; i <= n; ++i)
	{
	  He_nm2 = He_nm1;
	  He_nm1 = He_n;
	  He_n = x * He_nm1 - Tp(i - 1) * He_nm2;
	}

      return {n, x, He_n, He_nm1, He_nm2};
    }

  /**
   * Build a vector of the Gauss-Hermite integration rule abscissae and weights.
   */
  template<typename Tp>
    std::vector<emsr::QuadraturePoint<Tp>>
    hermite_zeros(unsigned int n, Tp proto = Tp{})
    {
      const auto s_eps = emsr::epsilon(proto);
      const unsigned int s_maxit = 1000u;
      const auto s_pim4 = Tp{0.7511255444649424828587030047762276930510L};
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;

      std::vector<emsr::QuadraturePoint<Tp>> pt(n);

      const auto m = n / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (n & 1)
	{
	  if (n < s_num_factorials<Tp>)
	    {
	      auto nm = n - 1;
	      auto nmfact = factorial<Tp>(nm);
	      auto mm = nm / 2;
	      auto mmfact = factorial<Tp>(mm);
	      auto Hnm1 = (mm & 1 ? Tp{-1} : Tp{1}) / mmfact;
	      pt[m].point = Tp{0};
	      pt[m].weight = s_sqrt_pi * std::pow(Tp{2}, Tp(n - 1))
				 / nmfact / Hnm1 / Hnm1 / n;
	    }
	  else
	    {
	      auto nm = n - 1;
	      auto nmfact = log_factorial<Tp>(nm);
	      auto mm = nm / 2;
	      auto mmfact = log_factorial<Tp>(mm);
	      pt[m].point = Tp{0};
	      pt[m].weight = s_sqrt_pi * std::pow(Tp{2}, Tp(n - 1))
				 *std::exp(-(nmfact - 2 * mmfact)) / n;
	    }
	}

      for (auto i = 1u; i <= m; ++i)
	{
	  Tp z;
	  Tp w = Tp{0};
	  if (i == 1)
	    z = std::sqrt(Tp(2 * n + 1))
		- 1.85575 * std::pow(Tp(2 * n + 1), -0.166667);
	  else if (i == 2)
	    z -= 1.14 * std::pow(Tp(n), 0.426) / z;
	  else if (i == 3)
	    z = 1.86 * z - 0.86 * pt[0].point;
	  else if (i == 4)
	    z = 1.91 * z - 0.91 * pt[1].point;
	  else
	    z = 2.0 * z - pt[i - 3].point;
	  for (auto its = 1u; its <= s_maxit; ++its)
	    {
	      auto H = s_pim4;
	      auto H1 = Tp{0};
	      for (auto k = 1u; k <= n; ++k)
		{
		  auto H2 = H1;
		  H1 = H;
		  H = z * std::sqrt(Tp{2} / k) * H1
		       - std::sqrt(Tp(k - 1) / Tp(k)) * H2;
		}
	      auto Hp = std::sqrt(Tp(2 * n)) * H1;
	      auto z1 = z;
	      z = z1 - H / Hp;
	      if (std::abs(z - z1) <= s_eps)
		{
		  w = Tp{2} / (Hp * Hp);
		  break;
		}
	      if (its > s_maxit)
		throw std::logic_error("hermite_zeros: Too many iterations");
	    }
	  pt[n - i].point = -z;
	  pt[n - i].weight = w;
	  pt[i - 1].point = z;
	  pt[i - 1].weight = w;
	}

      return pt;
    }

} // namespace detail
} // namespace emsr

#endif // SF_HERMITE_TCC
