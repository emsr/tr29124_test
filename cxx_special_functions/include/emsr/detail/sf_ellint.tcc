
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

/** @file bits/sf_ellint.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland.
//
// References:
// (1)  B. C. Carlson Numer. Math. 33, 1 (1979)
// (2)  B. C. Carlson, Special Functions of Applied Mathematics (1977)
// (3)  The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (4)  Numerical Recipes in C, 2nd ed, by W. H. Press, S. A. Teukolsky,
//      W. T. Vetterling, B. P. Flannery, Cambridge University Press
//      (1992), pp. 261-269
// (5)  Toshio Fukushima, Elliptic functions and elliptic integrals for
//      celestial mechanics and dynamical astronomy

#ifndef SF_ELLINT_TCC
#define SF_ELLINT_TCC 1

#include <stdexcept>
#include <complex>

#include <emsr/complex_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>

namespace emsr
{
namespace detail
{
  /**
   * @brief  Return the Carlson elliptic function
   * 	     @f$ R_C(x,y) = R_F(x,y,y) @f$ where @f$ R_F(x,y,z) @f$
   * 	     is the Carlson elliptic function of the first kind.
   *
   * The Carlson elliptic function is defined by:
   * @f[
   * 	 R_C(x,y) = \frac{1}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first argument.
   * @param  y  The second argument.
   * @return  The Carlson elliptic function.
   */
  template<typename Tp>
    Tp
    ellint_rc(Tp x, Tp y)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_min = emsr::lim_min(Real{});
      const auto s_eps = emsr::epsilon(Real{});
      const auto s_lolim = Real{5} * s_min;

      bool neg_x = false, neg_y = false;
      if constexpr (!emsr::is_complex_v<Tp>)
	{
	  if (std::real(x) < Real{0})
	    neg_x = true;
	  if (std::real(y) < Real{0})
	    neg_y = true;
	}

      if (std::isnan(x) || std::isnan(y))
	return s_NaN;
      else if (neg_x)
	throw std::domain_error("ellint_rc: argument less than zero");
      else if (std::abs(x) + std::abs(y) < s_lolim)
        throw std::domain_error("ellint_rc: arguments too small");
      else if (neg_y)
	{
	  if (std::abs(x) == Real{0})
	    return Tp{};
	  else
	    return std::sqrt(x / (x - y)) * ellint_rc(x - y, -y);
	}
      else
	{
	  auto xt = x;
	  auto yt = y;
	  auto A0 = (x + Real{2} * y) / Real{3};
	  auto Q = std::pow(Real{3} * s_eps, -Real{1} / Real{8})
		  * std::abs(A0 - x);
	  auto A = A0;
	  auto f = Real{1};

	  while (true)
	    {
	      auto lambda = Real{2} * std::sqrt(xt) * std::sqrt(yt)
			    + yt;
	      A = (A + lambda) / Real{4};
	      xt = (xt + lambda) / Real{4};
	      yt = (yt + lambda) / Real{4};
	      f *= Real{4};
	      if (Q < f * std::abs(A))
		{
		  auto s = (y - A0) / (f * A);
		  return (Real{1} + s * s * (Real{3} / Real{10}
			+ s * (Real{1} / Real{7}
			+ s * (Real{3} / Real{8}
			+ s * (Real{9} / Real{22}
			+ s * (Real{159} / Real{208}
			+ s * (Real{9} / Real{8}))))))) / std::sqrt(A);
		}
	    }

	  return Tp{0};
	}
    }

  /**
   * @brief  Return the Carlson elliptic function of the second kind
   * 	     @f$ R_D(x,y,z) = R_J(x,y,z,z) @f$ where
   * 	     @f$ R_J(x,y,z,p) @f$ is the Carlson elliptic function
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the second kind is defined by:
   * @f[
   * 	 R_D(x,y,z) = \frac{3}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{3/2}}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of two symmetric arguments.
   * @param  y  The second of two symmetric arguments.
   * @param  z  The third argument.
   * @return  The Carlson elliptic function of the second kind.
   */
  template<typename Tp>
    Tp
    ellint_rd(Tp x, Tp y, Tp z)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_min = emsr::lim_min(Real{});
      const auto s_eps = emsr::epsilon(Real{});
      const auto s_lolim = Real{5} * s_min;

      bool neg_arg = false;
      if constexpr (!emsr::is_complex_v<Tp>)
	if (std::real(x) < Real{0}
	 || std::real(y) < Real{0}
	 || std::real(z) < Real{0})
	  neg_arg = true;

      if (std::isnan(x) || std::isnan(y) || std::isnan(z))
	return s_NaN;
      if (neg_arg)
	throw std::domain_error("ellint_rd: argument less than zero");
      else if (std::abs(x) + std::abs(y) < s_lolim
	    || std::abs(z) < s_lolim)
	throw std::domain_error("ellint_rd: arguments too small");
      else
	{
	  auto xt = x;
	  auto yt = y;
	  auto zt = z;
	  auto A0 = (x + y + Real{3} * z) / Real{5};
	  auto Q = std::pow(s_eps / Real{4}, -Real{1} / Real{6})
		  * std::max(std::abs(A0 - z),
			     std::max(std::abs(A0 - x),
				      std::abs(A0 - y)));
	  auto A = A0;
	  auto f = Real{1};
	  auto sum = Tp{0};

	  while (true)
	    {
	      auto lambda = std::sqrt(xt) * std::sqrt(yt)
			    + std::sqrt(yt) * std::sqrt(zt)
			    + std::sqrt(zt) * std::sqrt(xt);
	      sum += Real{1} / f / std::sqrt(zt) / (zt + lambda);
	      A = (A + lambda) / Real{4};
	      xt = (xt + lambda) / Real{4};
	      yt = (yt + lambda) / Real{4};
	      zt = (zt + lambda) / Real{4};
	      f *= Real{4};
	      if (Q < f * std::abs(A))
		{
		  auto Xi = (A0 - x) / (f * A);
		  auto Yi = (A0 - y) / (f * A);
		  auto Zi = -(Xi + Yi) / Real{3};
		  auto ZZ = Zi * Zi;
		  auto XY = Xi * Yi;
		  auto E2 = XY - Real{6} * ZZ;
		  auto E3 = (Real{3} * XY - Real{8} * ZZ) * Zi;
		  auto E4 = Real{3} * (XY - ZZ) * ZZ;
		  auto E5 = XY * Zi * ZZ;
		  return (Real{1}
			- Real{3} * E2 / Real{14}
			+ E3 / Real{6}
			+ Real{9} * E2 * E2 / Real{88}
			- Real{3} * E4 / Real{22}
			- Real{9} * E2 * E3 / Real{52}
			+ Real{3} * E5 / Real{26}) / f / A / std::sqrt(A)
			+ Real{3} * sum;
		}
	    }

	  return Tp{0};
	}
    }

  template<typename Tp>
    Tp
    comp_ellint_rf(Tp x, Tp y)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_tolfact = Real{2.7L} * emsr::sqrt_eps(Real{});

      if (std::isnan(x) || std::isnan(y))
	return s_NaN;
      else
	{
	  x = std::sqrt(x);
	  y = std::sqrt(y);
	  while (true)
	    {
	      auto xt = x;
	      x = (x + y) / Real{2};
	      y = std::sqrt(xt) * std::sqrt(y);
	      if (std::abs(x - y) < s_tolfact * std::abs(x))
		return s_pi / (x + y);
	    }
	}
    }

  /**
   * @brief Return the Carlson elliptic function @f$ R_F(x,y,z) @f$
   * 	    of the first kind.
   *
   * The Carlson elliptic function of the first kind is defined by:
   * @f[
   * 	 R_F(x,y,z) = \frac{1}{2} \int_0^\infty
   * 		   \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}}
   * @f]
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   * @return  The Carlson elliptic function of the first kind.
   */
  template<typename Tp>
    Tp
    ellint_rf(Tp x, Tp y, Tp z)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_min = emsr::lim_min(Real{});
      const auto s_eps = emsr::epsilon(Real{});
      const auto s_lolim = Real(5) * s_min;

      bool neg_arg = false;
      if constexpr (!emsr::is_complex_v<Tp>)
	if (std::real(x) < Real{0}
	 || std::real(y) < Real{0}
	 || std::real(z) < Real{0})
	  neg_arg = true;

      if (std::isnan(x) || std::isnan(y) || std::isnan(z))
	return s_NaN;
      else if (neg_arg)
        throw std::domain_error("ellint_rf: argument less than zero");
      else if (std::abs(x) + std::abs(y) < s_lolim
	    || std::abs(x) + std::abs(z) < s_lolim
	    || std::abs(y) + std::abs(z) < s_lolim)
        throw std::domain_error("ellint_rf: argument too small");

      if (std::abs(z) < s_eps)
        return comp_ellint_rf(x, y);
      else if (std::abs(z - y) < s_eps)
	return ellint_rc(x, y);
      else
	{
	  auto xt = x;
	  auto yt = y;
	  auto zt = z;
	  auto A0 = (x + y + z) / Real{3};
	  auto Q = std::pow(Real{3} * s_eps, -Real{1} / Real{6})
		  * std::max(std::abs(A0 - z),
			     std::max(std::abs(A0 - x),
				      std::abs(A0 - y)));
	  auto A = A0;
	  auto f = Real{1};

	  while (true)
	    {
	      auto lambda = std::sqrt(xt) * std::sqrt(yt)
			    + std::sqrt(yt) * std::sqrt(zt)
			    + std::sqrt(zt) * std::sqrt(xt);
	      A = (A + lambda) / Real{4};
	      xt = (xt + lambda) / Real{4};
	      yt = (yt + lambda) / Real{4};
	      zt = (zt + lambda) / Real{4};
	      f *= Real{4};
	      if (Q < f * std::abs(A))
		{
		  auto Xi = (A0 - x) / (f * A);
		  auto Yi = (A0 - y) / (f * A);
		  auto Zi = -(Xi + Yi);
		  auto E2 = Xi * Yi - Zi * Zi;
		  auto E3 = Xi * Yi * Zi;
		  return (Real{1}
			- E2 / Real{10}
			+ E3 / Real{14}
			+ E2 * E2 / Real{24}
			- Real{3} * E2 * E3 / Real{44}) / std::sqrt(A);
		}
	    }

	  return Tp{0};
	}
    }

  template<typename Tp>
    Tp
    comp_ellint_rg(Tp x, Tp y)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_pi = emsr::pi_v<Real>;
      const auto s_tolfact = Real{2.7L} * emsr::sqrt_eps(Real{});

      if (std::isnan(x) || std::isnan(y))
	return s_NaN;
      else if (x == Tp{0} && y == Tp{})
	return Tp{};
      else if (x == Tp{0})
	return std::sqrt(y) / Real{2};
      else if (y == Tp{0})
	return std::sqrt(x) / Real{2};
      else
	{
	  auto xt = std::sqrt(x);
	  auto yt = std::sqrt(y);
	  const auto A = (xt + yt) / Real{2};
	  auto sum = Tp{};
	  auto sf = Real{1} / Real{2};
	  while (true)
	    {
	      auto xtt = xt;
	      xt = (xt + yt) / Real{2};
	      yt = std::sqrt(xtt) * std::sqrt(yt);
	      auto del = xt - yt;
	      if (std::abs(del) < s_tolfact * std::abs(xt))
		return (A * A - sum) * s_pi / (xt + yt) / Real{2};
	      sum += sf * del * del;
	      sf *= Real{2};
	    }
	}
    }

  /**
   * @brief  Return the symmetric Carlson elliptic function of the second kind
   * 	     @f$ R_G(x,y,z) @f$.
   *
   * The Carlson symmetric elliptic function of the second kind is defined by:
   * @f[
   * 	 R_G(x,y,z) = \frac{1}{4} \int_0^\infty
   * 		   dt t [(t + x)(t + y)(t + z)]^{-1/2}
   * 		   (\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z})
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   * @return  The Carlson symmetric elliptic function of the second kind.
   */

  template<typename Tp>
    Tp
    ellint_rg(Tp x, Tp y, Tp z)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});

      if (std::isnan(x) || std::isnan(y) || std::isnan(z))
	return s_NaN;
      else if (z == Tp{0})
	return comp_ellint_rg(x, y);
      else if (x == Tp{0})
	return comp_ellint_rg(y, z);
      else if (y == Tp{0})
	return comp_ellint_rg(z, x);
      else
	//return (z * ellint_rf(x, y, z)
	//     - (x - z) * (y - z) * ellint_rd(x, y, z) / Real{3}
	//     + (std::sqrt(x) * std::sqrt(y) / std::sqrt(z))) / Real{2};
	// There is a symmetric version that is less subject to cancellation loss
	// when the arguments are real:
	return (x * (y + z) * ellint_rd(y, z, x)
	      + y * (z + x) * ellint_rd(z, x, y)
	      + z * (x + y) * ellint_rd(x, y, z)) / Real{6};
    }

  /**
   * @brief  Return the Carlson elliptic function @f$ R_J(x,y,z,p) @f$
   * 	     of the third kind.
   *
   * The Carlson elliptic function of the third kind is defined by:
   * @f[
   * 	 R_J(x,y,z,p) = \frac{3}{2} \int_0^\infty
   * 	 \frac{dt}{(t + x)^{1/2}(t + y)^{1/2}(t + z)^{1/2}(t + p)}
   * @f]
   *
   * Based on Carlson's algorithms:
   * -  B. C. Carlson Numer. Math. 33, 1 (1979)
   * -  B. C. Carlson, Special Functions of Applied Mathematics (1977)
   * -  Numerical Recipes in C, 2nd ed, pp. 261-269,
   * 	by Press, Teukolsky, Vetterling, Flannery (1992)
   *
   * @param  x  The first of three symmetric arguments.
   * @param  y  The second of three symmetric arguments.
   * @param  z  The third of three symmetric arguments.
   * @param  p  The fourth argument.
   * @return  The Carlson elliptic function of the fourth kind.
   */
  template<typename Tp>
    Tp
    ellint_rj(Tp x, Tp y, Tp z, Tp p)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});
      const auto s_min = emsr::lim_min(Real{});
      const auto s_eps = emsr::epsilon(Real{});
      const auto s_lolim = Real(5) * s_min;

      bool neg_arg = false;
      if constexpr (!emsr::is_complex_v<Tp>)
	if (std::real(x) < Real{0}
	 || std::real(y) < Real{0}
	 || std::real(z) < Real{0})
	  neg_arg = true;

      if (std::isnan(x) || std::isnan(y) || std::isnan(z) || std::isnan(p))
	return s_NaN;
      else if (neg_arg)
        throw std::domain_error("ellint_rj: argument less than zero");
      else if (std::abs(x) + std::abs(y) < s_lolim
	    || std::abs(y) + std::abs(z) < s_lolim
	    || std::abs(z) + std::abs(x) < s_lolim
	    || std::abs(p) < s_lolim)
        throw std::domain_error("ellint_rj: argument too small");
      else if (std::abs(p - z) < s_eps)
	return ellint_rd(x, y, z);
      else
	{
	  auto xt = x;
	  auto yt = y;
	  auto zt = z;
	  auto pt = p;
	  auto A0 = (x + y + z + Real{2} * p) / Real{5};
	  auto delta = (p - x) * (p - y) * (p - z);
	  auto Q = std::pow(s_eps / Real{4}, -Real{1} / Real{6})
		  * std::max(std::abs(A0 - z),
		      std::max(std::abs(A0 - x),
			std::max(std::abs(A0 - y), std::abs(A0 - p))));
	  auto A = A0;
	  auto f = Real{1};
	  auto fe = Real{1};
	  auto sum = Tp{0};

	  while (true)
	    {
	      auto xroot = std::sqrt(xt);
	      auto yroot = std::sqrt(yt);
	      auto zroot = std::sqrt(zt);
	      auto proot = std::sqrt(pt);
	      auto lambda = xroot * yroot
			    + yroot * zroot
			    + zroot * xroot;
	      A = (A + lambda) / Real{4};
	      xt = (xt + lambda) / Real{4};
	      yt = (yt + lambda) / Real{4};
	      zt = (zt + lambda) / Real{4};
	      pt = (pt + lambda) / Real{4};
	      auto d = (proot + xroot)
		       * (proot + yroot)
		       * (proot + zroot);
	      auto E = delta / (fe * d * d);
	      sum += ellint_rc(Tp{1}, Tp{1} + E) / (f * d);
	      f *= Real{4};
	      fe *= Real{64};
	      if (Q < f * std::abs(A))
		{
		  auto Xi = (A0 - x) / (f * A);
		  auto Yi = (A0 - y) / (f * A);
		  auto Zi = (A0 - z) / (f * A);
		  auto XYZ = Xi * Yi * Zi;
		  auto Pi = -(Xi + Yi + Zi) / Real{2};
		  auto PP = Pi * Pi;
		  auto PPP = PP * Pi;
		  auto E2 = Xi * Yi
			   + Yi * Zi
			   + Zi * Xi
			   - Real{3} * PP;
		  auto E3 = XYZ + Real{2} * E2 * Pi + Tp{4} * PPP;
		  auto E4 = Pi
			   * (Real{2} * XYZ + E2 * Pi + Real{3} * PPP);
		  auto E5 = XYZ * PP;
		  return (Real{1} - Real{3} * E2 / Real{14}
			+ E3 / Real{6}
			+ Real{9} * E2 * E2 / Real{88}
			- Real{3} * E4 / Real{22}
			- Real{9} * E2 * E3 / Real{52}
			+ Real{3} * E5 / Real{26}) / f / A / std::sqrt(A)
			+ Real{6} * sum;
		}
	    }

	  return Tp{};
	}
    }

  /**
   * @brief  Return the complete elliptic integral of the first kind
   * @f$ K(k) @f$ using the Carlson formulation.
   *
   * The complete elliptic integral of the first kind is defined as
   * @f[
   *   K(k) = F(k,\pi/2) = \int_0^{\pi/2}\frac{d\theta}
   * 					     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   * where @f$ F(k,\phi) @f$ is the incomplete elliptic integral of the
   * first kind.
   *
   * @param  k  The modulus of the complete elliptic function.
   * @return  The complete elliptic function of the first kind.
   */
  template<typename Tp>
    Tp
    comp_ellint_1(Tp k)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);

      if (std::isnan(k))
	return s_NaN;
      else if (std::abs(k) == Real{1})
	return s_NaN;
      else
	return comp_ellint_rf(Tp{1} - k * k, Tp{1});
    }

  /**
   * @brief  Return the incomplete elliptic integral of the first kind
   * @f$ F(k,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the first kind is defined as
   * @f[
   *   F(k,\phi) = \int_0^{\phi}\frac{d\theta}
   * 				     {\sqrt{1 - k^2 sin^2\theta}}
   * @f]
   *
   * @param  k  The elliptic modulus.
   * @param  phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the first kind.
   */
  template<typename Tp>
    Tp
    ellint_1(Tp k, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);
      const auto s_pi = emsr::pi_v<Real>;

      if (std::isnan(k) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Real{1})
	throw std::domain_error("ellint_1: bad argument");
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(std::real(phi) / s_pi + Real{0.5L});
	  const auto phi_red = phi - n * s_pi;

	  const auto s = std::sin(phi_red);
	  const auto c = std::cos(phi_red);

	  const auto F = s
			 * ellint_rf(c * c,
				       Tp{1} - k * k * s * s,
				       Tp{1});

	  if (n == 0)
	    return F;
	  else
	    return F + Tp{2} * n * comp_ellint_1(k);
	}
    }

  /**
   * @brief  Return the complete elliptic integral of the second kind
   * 	     @f$ E(k) @f$ using the Carlson formulation.
   *
   * The complete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\pi/2) = \int_0^{\pi/2}\sqrt{1 - k^2 sin^2\theta}
   * @f]
   *
   * @param  k  The modulus of the complete elliptic function.
   * @return  The complete elliptic function of the second kind.
   */
  template<typename Tp>
    Tp
    comp_ellint_2(Tp k)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);

      if (std::isnan(k))
	return s_NaN;
      else if (std::abs(k) == Real{1})
	return Tp{1};
      else if (std::abs(k) > Real{1})
	throw std::domain_error("comp_ellint_2: bad argument");
      else
	{
	  const auto kk = k * k;

	  return ellint_rf(Tp{0}, Tp{1} - kk, Tp{1})
	       - kk * ellint_rd(Tp{0}, Tp{1} - kk, Tp{1})
	       / Tp{3};
	}
    }

  /**
   * @brief  Return the incomplete elliptic integral of the second kind
   * @f$ E(k,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the second kind is defined as
   * @f[
   *   E(k,\phi) = \int_0^{\phi} \sqrt{1 - k^2 sin^2\theta}
   * @f]
   *
   * @param  k  The elliptic modulus.
   * @param  phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the second kind.
   */
  template<typename Tp>
    Tp
    ellint_2(Tp k, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);
      const auto s_pi = emsr::pi_v<Real>;

      if (std::isnan(k) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Real{1})
	throw std::domain_error("ellint_2: bad argument");
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(std::real(phi) / s_pi + Real{0.5L});
	  const auto phi_red = phi - n * s_pi;

	  const auto kk = k * k;
	  const auto s = std::sin(phi_red);
	  const auto ss = s * s;
	  const auto sss = ss * s;
	  const auto c = std::cos(phi_red);
	  const auto cc = c * c;

	  const auto E = s
		        * ellint_rf(cc, Tp{1} - kk * ss, Tp{1})
		        - kk * sss
		        * ellint_rd(cc, Tp{1} - kk * ss, Tp{1})
		        / Tp{3};

	  if (n == 0)
	    return E;
	  else
	    return E + Tp{2} * n * comp_ellint_2(k);
	}
    }

  /**
   * @brief Return the complete elliptic integral of the third kind
   * @f$ \Pi(k,\nu) = \Pi(k,\nu,\pi/2) @f$ using the
   * Carlson formulation.
   *
   * The complete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu) = \int_0^{\pi/2}
   * 		     \frac{d\theta}
   * 		   {(1 - \nu \sin^2\theta)\sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   *
   * @param  k  The elliptic modulus.
   * @param  nu  The characteristic.
   * @return  The complete elliptic function of the third kind.
   */
  template<typename Tp>
    Tp
    comp_ellint_3(Tp k, Tp nu)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);

      if (std::isnan(k) || std::isnan(nu))
	return s_NaN;
      else if (nu == Tp{1})
	return emsr::infinity(k);
      else if (std::abs(k) > Real{1})
	throw std::domain_error("comp_ellint_3: bad argument");
      else
	{
	  const auto mc = Tp{1} - k * k;

	  return ellint_rf(Tp{0}, mc, Tp{1})
	     + nu * ellint_rj(Tp{0}, mc, Tp{1}, Tp{1} - nu) / Tp{3};
	}
    }

  /**
   * @brief Return the incomplete elliptic integral of the third kind
   * @f$ \Pi(k,\nu,\phi) @f$ using the Carlson formulation.
   *
   * The incomplete elliptic integral of the third kind is defined as
   * @f[
   *   \Pi(k,\nu,\phi) = \int_0^{\phi}
   * 			 \frac{d\theta}
   * 			      {(1 - \nu \sin^2\theta)
   * 			       \sqrt{1 - k^2 \sin^2\theta}}
   * @f]
   *
   * @param  k  The elliptic modulus.
   * @param  nu  The characteristic.
   * @param  phi  The integral limit argument of the elliptic function.
   * @return  The elliptic function of the third kind.
   */
  template<typename Tp>
    Tp
    ellint_3(Tp k, Tp nu, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);
      const auto s_pi = emsr::pi_v<Real>;

      if (std::isnan(k) || std::isnan(nu) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Real{1})
	throw std::domain_error("ellint_3: bad argument");
      else
	{
	  // Reduce phi to -pi/2 < phi < +pi/2.
	  const int n = std::floor(std::real(phi) / s_pi + Real{0.5L});
	  const auto phi_red = phi - n * s_pi;

	  const auto kk = k * k;
	  const auto s = std::sin(phi_red);
	  const auto ss = s * s;
	  const auto sss = ss * s;
	  const auto c = std::cos(phi_red);
	  const auto cc = c * c;

	  const auto Pi = s
			 * ellint_rf(cc, Tp{1} - kk * ss, Tp{1})
			 + nu * sss
			 * ellint_rj(cc, Tp{1} - kk * ss, Tp{1},
				       Tp{1} - nu * ss) / Tp{3};

	  if (n == 0)
	    return Pi;
	  else
	    return Pi + Tp{2} * n * comp_ellint_3(k, nu);
	}
    }

  /**
   * Return the Legendre elliptic integral D.
   */
  template<typename Tp>
    Tp
    ellint_d(Tp k, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);

      if (std::isnan(k) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Real{1})
	throw std::domain_error("ellint_d: bad argument");
      else
	{
	  auto sinphi = std::sin(phi);
	  auto sinphi2 = sinphi * sinphi;
	  auto k2 = k * k;
	  auto arg1 = Tp{1} - sinphi2;
	  auto arg2 = Tp{1} - k2 * sinphi2;
	  return sinphi * sinphi2 * ellint_rd(arg1, arg2, Tp{1})
		 / Tp{3};
	}
    }

  /**
   * Return the complete Legendre elliptic integral D.
   */
  template<typename Tp>
    Tp
    comp_ellint_d(Tp k)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);

      if (std::isnan(k))
	return s_NaN;
      else
	return ellint_rd(Tp{0}, Tp{1} - k * k, Tp{1}) / Real{3};
    }

  /**
   * Return the Bulirsch elliptic integrals of the first kind.
   */
  template<typename Tp>
    Tp
    ellint_el1(Tp x, Tp k_c)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});

      if (std::isnan(x) || std::isnan(k_c))
	return s_NaN;
      else
	{
	  auto x2 = x * x;
	  auto k2_c = k_c * k_c;
	  auto arg2 = Tp{1} + k2_c * x2;
	  return x * ellint_rf(Tp{1}, arg2, Tp{1} + x2);
	}
    }

  /**
   * Return the Bulirsch elliptic integrals of the second kind.
   */
  template<typename Tp>
    Tp
    ellint_el2(Tp x, Tp k_c, Tp a, Tp b)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});

      if (std::isnan(x) || std::isnan(k_c) || std::isnan(a) || std::isnan(b))
	return s_NaN;
      else
	{
	  auto x2 = x * x;
	  auto x3 = x * x;
	  auto k2_c = k_c * k_c;
	  auto arg2 = Tp{1} + k2_c * x2;
	  auto arg3 = Tp{1} + x2;
	  return a * x * ellint_rf(Tp{1}, arg2, arg3)
               + (b - a) * x3
		 * ellint_rd(Tp{1}, arg2, arg3) / Tp{3};
	}
    }

  /**
   * Return the Bulirsch elliptic integrals of the third kind.
   */
  template<typename Tp>
    Tp
    ellint_el3(Tp x, Tp k_c, Tp p)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(Real{});

      if (std::isnan(x) || std::isnan(k_c) || std::isnan(p))
	return s_NaN;
      else
	{
	  auto x2 = x * x;
	  auto x3 = x * x;
	  auto k2_c = k_c * k_c;
	  auto arg2 = Tp{1} + k2_c * x2;
	  auto arg3 = Tp{1} + x2;
	  auto arg4 = Tp{1} + p * x2;
	  return x * ellint_rf(Real{1}, arg2, arg3)
	       + (Tp{1} - p) * x3
		 * ellint_rj(Tp{1}, arg2, arg3, arg4) / Tp{3};
	}
    }

  /**
   * Return the Bulirsch complete elliptic integrals.
   */
  template<typename Tp>
    Tp
    ellint_cel(Tp k_c, Tp p, Tp a, Tp b)
    {
      const auto s_NaN = emsr::quiet_NaN(k_c);

      if (std::isnan(k_c) || std::isnan(p) || std::isnan(a) || std::isnan(b))
	return s_NaN;
      else
	{
	  auto k2_c = k_c * k_c;
	  return a * ellint_rf(Tp{0}, k2_c, Tp{1})
	       + (b - p * a)
		 * ellint_rj(Tp{0}, k2_c, Tp{1}, p) / Tp{3};
	}
    }

  /**
   * Return the Jacobi zeta function.
   */
  template<typename Tp>
    Tp
    jacobi_zeta(Tp k, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_eps = emsr::epsilon(k);

      if (std::isnan(k) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Real{1})
	throw std::domain_error("jacobi_zeta: bad argument");
      else if (std::abs(k) < s_eps)
	return Tp{0};
      else if (std::abs(k - Real{1}) < s_eps)
	return std::sin(phi);
      else if (std::real(phi) < Real{0})
	return -jacobi_zeta(k, -phi);
      else if (std::abs(phi) < s_eps || std::abs(phi - s_pi_2) < s_eps)
	return Tp{0};
      else
	{
	  auto mc = k * k;
	  auto cosphi = std::cos(phi);
	  auto sinphi = std::sin(phi);
	  auto m = Tp{1} - mc;
	  auto arg4 = Tp{1} - mc  * sinphi * sinphi;
	  return mc * cosphi * sinphi * std::sqrt(arg4)
	       * ellint_rj(Tp{0}, m, Real{1}, arg4)
	       / (Tp{3} * comp_ellint_1(k));
	}
    }

  /**
   * Return the Heuman lambda function.
   */
  template<typename Tp>
    Tp
    heuman_lambda(Tp k, Tp phi)
    {
      using Real = emsr::num_traits_t<Tp>;
      const auto s_NaN = emsr::quiet_NaN(k);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
      const auto s_eps = emsr::epsilon(k);

      if (std::isnan(k) || std::isnan(phi))
	return s_NaN;
      else if (std::abs(k) > Tp{1})
	throw std::domain_error("heuman_lambda: bad argument");
      else if (std::abs(std::abs(k) - Tp{1}) < s_eps)
	return phi / s_pi_2;
      else if (std::abs(k) < s_eps)
	return std::sin(phi);
      else if (std::real(phi) < Real{0})
	return -heuman_lambda(k, -phi);
      else if (std::abs(phi - s_pi_2) < Tp{5} * s_eps)
	return Tp{1};
      else if (std::abs(phi) < Tp{5} * s_eps)
	return Tp{0};
      else if (std::abs(phi) < s_pi_2)
	{
	  auto mc = k * k;
	  auto m = Tp{1} - mc;
	  auto cosphi = std::cos(phi);
	  auto sinphi = std::sin(phi);
	  auto Delta2 = Tp{1} - m * sinphi * sinphi;
	  if (std::abs(Delta2) < Real{0})
	    Delta2 = Tp{0};
	  auto fact = Tp{2} * m * cosphi * sinphi
		      / (s_pi * std::sqrt(Delta2));
	  auto arg4 = Tp{1} - mc / Delta2;
	  auto fact2 = mc / (Tp{3} * Delta2);

	  return fact * (ellint_rf(Tp{0}, m, Tp{1})
			+ fact2 * ellint_rj(Tp{0}, m, Tp{1}, arg4));
	}
      else
	{
	  auto kc = std::sqrt(Tp{1} - k * k);
	  return ellint_1(kc, phi) / comp_ellint_1(kc)
	       + comp_ellint_1(k) * jacobi_zeta(kc, phi) / s_pi_2;
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_ELLINT_TCC

