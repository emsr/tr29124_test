
// Copyright (C) 2006-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/specfun_state.h
 * This is an internal header file, included by other library headers.
 * Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SPECFUN_STATE_H
#define SPECFUN_STATE_H 1

#include <cmath>
#include <vector>
#include <limits>

namespace emsr
{

  /** 
   * \brief Enumeration for differing types of Gauss quadrature.
   *        The gauss_quad_type is used to determine the boundary condition
   *        modifications applied to orthogonal polynomials for quadrature
   *        rules.
   */
  enum gauss_quad_type
  {
    Gauss,             ///< Gauss quadrature
    Gauss_Lobatto,     ///< Gauss-Lobatto quadrature
    Gauss_Radau_lower, ///< Gauss-Radau quadrature including the node -1
    Gauss_Radau_upper  ///< Gauss-Radau quadrature including the node +1
  };

  /**
   * A type describing the state of a Hermite polynomial.
   */
  template<typename Tp>
    struct hermite_t
    {
      unsigned int n;
      Tp x;
      Tp H_n;
      Tp H_nm1;
      Tp H_nm2;

      constexpr Tp
      deriv() const noexcept
      { return Tp(2 * n) * H_nm1; }

      constexpr Tp
      deriv2() const noexcept
      { return Tp(4 * n * (n - 1)) * H_nm2; }
    };

  /**
   * A type describing the state of a probabilists Hermite polynomial.
   */
  template<typename Tp>
    struct hermite_he_t
    {
      unsigned int n;
      Tp x;
      Tp He_n;
      Tp He_nm1;
      Tp He_nm2;

      constexpr Tp
      deriv() const noexcept
      { return Tp(n) * He_nm1; }

      constexpr Tp
      deriv2() const noexcept
      { return Tp(n * (n - 1)) * He_nm2; }
    };

  /**
   * A type describing the state of a Legendre polynomial.
   *
   * The method lobatto() will return the Lobatto polynomial:
   * @f[
   *   Lo_l(x) = (1 - x^2)P'_l(x) = l \left[ P_{l - 1}(x) - x P_l(x) \right]
   * @f]
   */
  template<typename Tp>
    struct legendre_p_t
    {
      unsigned int l;
      Tp x;
      Tp P_l;   /// P_l(x)
      Tp P_lm1; /// P_{l-1}(x)
      Tp P_lm2; /// P_{l-2}(x)

      // Return the Lobatto polynomial.
      constexpr Tp
      lobatto() const noexcept
      { return l * (P_l - x * P_lm1); }

      constexpr Tp
      deriv() const noexcept
      {
	if (std::abs(x) == Tp{1})
	  {
	    const auto sgn = x == Tp{+1}
			     ? Tp{+1}
			     : (l % 2 == 0 ? Tp{-1} : Tp{+1});
	    return sgn * Tp(l) * Tp(l + 1) / Tp{2};
	  }
	else
	  return l * (x * P_l - P_lm1)
	       / ((Tp{1} - x) * (Tp{1} + x));
      }
    };

  /**
   * A type describing the state of an associated Legendre function.
   */
  template<typename Tp>
    struct assoc_legendre_p_t
    {
      unsigned int l;
      unsigned int m;
      Tp x;
      Tp P_lm;   /// P_l^{(m)}(x)
      Tp P_lm1m; /// P_{l-1}^{(m)}(x)
      Tp P_lm2m; /// P_{l-2}^{(m)}(x)
      Tp phase = 1; // -1 For Condon-Shortley.

      constexpr Tp
      deriv() const noexcept
      {
	if (std::abs(x) == Tp{1})
	  {
	    const auto sgn = x == Tp{+1}
			     ? Tp{+1}
			     : (l % 2 == 0 ? Tp{-1} : Tp{+1});
	    if (m == 0)
	      return sgn *  Tp(l) * Tp(l + 1) / Tp{2};
	    else if (m == 1)
	      {
		const auto sgn = x == Tp{+1}
				 ? Tp{+1}
				 : (l % 2 == 0 ? Tp{+1} : Tp{-1});
		return -phase * sgn * std::numeric_limits<Tp>::infinity();
	      }
	    else if (m == 2)
	      return -sgn * Tp(l + 2) * Tp(l + 1) / Tp{2}
		   * Tp(l) * Tp(int(l) - 1) / Tp{2};
	    else
	      return Tp{0};
	  }
	else
	  return -phase * ((l + m) * P_lm1m - l * x * P_lm)
	       / ((Tp{1} - x) * (Tp{1} + x));
      }
    };

  /**
   * A type describing the state of a Legendre function of the second kind.
   */
  template<typename Tp>
    struct legendre_q_t
    {
      unsigned int l;
      Tp x;
      Tp Q_l;   /// Q_l(x)
      Tp Q_lm1; /// Q_{l-1}(x)
      Tp Q_lm2; /// Q_{l-2}(x)

      constexpr Tp
      deriv() const noexcept
      {
	if (std::abs(x) == Tp{1})
	  return Tp(l % 2 == 1 ? -1 : +1)
		* std::numeric_limits<Tp>::infinity();
	else
	  return Tp(l) * (x * Q_l - Q_lm1)
	       / ((Tp{1} - x) * (Tp{1} + x));
      }
    };

  /**
   * A type describing the state of an associated Legendre function
   * of the second kind.
   */
  template<typename Tp>
    struct assoc_legendre_q_t
    {
      unsigned int l; /// degree
      unsigned int m; /// order
      Tp x; /// argument
      Tp Q_lm;   /// Q_l^{(m)}(x)
      Tp Q_lmm1; /// Q_l^{(m-1)}(x)
      Tp Q_lmm2; /// Q_l^{(m-2)}(x)
      Tp phase = 1; // -1 For Condon-Shortley.

      constexpr Tp
      deriv() const noexcept
      {
	if (std::abs(x) == 1)
	  return Tp(l % 2 == 1 ? -1 : +1)
		* std::numeric_limits<Tp>::infinity();
	else
	  {
	    const auto fact = (Tp{1} - x) * (Tp{1} + x);
	    const auto root = std::sqrt(Tp{1} - x)
			      * std::sqrt(Tp{1} + x);
	    return Tp(m) * x * Q_lm / fact
		 + Tp(l + m) * Tp(l - m + 1) * Q_lmm1 / root;
	  }
      }
    };

  /**
   * A type describing the state of a Laguerre polynomial.
   */
  template<typename Tpa, typename Tp>
    struct laguerre_t
    {
      unsigned int n;
      Tpa alpha1;
      Tp x;
      Tp L_n;
      Tp L_nm1;
      Tp L_nm2;

      constexpr Tp
      deriv() const noexcept
      { return (Tp(n) * L_nm1 - Tp(n + alpha1) * L_nm2) / x; }
    };

  /**
   * A type describing the state of a Jacobi polynomial.
   */
  template<typename Tp>
    struct jacobi_t
    {
      unsigned int n;
      Tp alpha1;
      Tp beta1;
      Tp x;
      Tp P_n;
      Tp P_nm1;
      Tp P_nm2;

      constexpr Tp
      deriv() const noexcept
      {
	auto apbp2k = alpha1 + beta1 + Tp(2 * n);
	return (n * (alpha1 - beta1 - apbp2k * x) * P_nm1
		   + Tp{2} * (n + alpha1) * (n + beta1) * P_nm2)
		/ (apbp2k * (Tp{1} - x * x));
      }
    };

  /**
   * A type describing the state of a Gegenbauer polynomial.
   */
  template<typename Tp>
    struct gegenbauer_t
    {
      unsigned int n;
      Tp lambda;
      Tp x;
      Tp C_n;
      Tp C_nm1;
      Tp C_nm2;

      constexpr Tp
      deriv() const noexcept
      {
	auto apbp2k = Tp{2} * lambda + Tp(2 * n);
	return (n * (-apbp2k * x) * C_nm1
		   + Tp{2} * (n + lambda) * (n + lambda) * C_nm2)
		/ (apbp2k * (Tp{1} - x * x));
      }
    };

  /**
   * A type describing the state of a Chebyshev polynomial of the first kind.
   */
  template<typename Tp>
    struct chebyshev_t_t
    {
      unsigned int n;
      Tp x;
      Tp T_n;
      Tp T_nm1;
      Tp T_nm2;

      constexpr Tp
      deriv() const noexcept
      { return Tp(n) * (T_nm1 - x * T_n) / (Tp{1} - x * x); }

      constexpr Tp
      deriv2() const noexcept
      {
	const auto xx = x * x;
	const auto num = Tp{1} - xx;
	return Tp(n)
	     * (x * T_nm1 + (Tp(n - 1) * xx - Tp(n)) * T_n)
	     / num / num;
      }
    };

  /**
   * A type describing the state of a Chebyshev polynomial of the second kind.
   */
  template<typename Tp>
    struct chebyshev_u_t
    {
      unsigned int n;
      Tp x;
      Tp U_n;
      Tp U_nm1;
      Tp U_nm2;

      constexpr Tp
      deriv() const noexcept
      {
	return (Tp(n + 1) * U_nm1 - Tp(n) * x * U_n)
		/ (Tp{1} - x * x);
      }
    };

  /**
   * A type describing the state of a Chebyshev polynomial of the third kind.
   */
  template<typename Tp>
    struct chebyshev_v_t
    {
      unsigned int n;
      Tp x;
      Tp V_n;
      Tp V_nm1;
      Tp V_nm2;

      constexpr Tp
      deriv() const noexcept
      {
	auto apbp2k = Tp(2 * n);
	return (n * (Tp{1} - apbp2k * x) * V_nm1
		   + Tp(2 * (n + 0.5L) * (n + -0.5L)) * V_nm2)
		/ (apbp2k * (Tp{1} - x * x));
      }
    };

  /**
   * A type describing the state of a Chebyshev polynomial of the fourth kind.
   */
  template<typename Tp>
    struct chebyshev_w_t
    {
      unsigned int n;
      Tp x;
      Tp W_n;
      Tp W_nm1;
      Tp W_nm2;

      constexpr Tp
      deriv() const noexcept
      {
	auto apbp2k = Tp(2 * n);
	return (n * (Tp{-1} - apbp2k * x) * W_nm1
		   + Tp(2 * (n - 0.5L) * (n + 0.5L)) * W_nm2)
		/ (apbp2k * (Tp{1} - x * x));
      }
    };

  template<typename Tx, typename Tp>
    struct airy_t
    {
      /// The argument of the Airy fuctions.
      Tx x_arg;

      /// The value of the Airy function Ai.
      Tp Ai_value;

      /// The derivative of the Airy function Ai.
      Tp Ai_deriv;

      /// The value of the Airy function Bi.
      Tp Bi_value;

      /// The derivative of the Airy function Bi.
      Tp Bi_deriv;

      /// Return the Wronskian of this Airy function state.
      constexpr Tp
      Wronskian() const noexcept
      { return Ai_value * Bi_deriv - Bi_value * Ai_deriv; }
    };

  /**
   * Tp pretty much has to be complex.
   */
  template<typename Tx, typename Tp>
    struct fock_airy_t
    {
      /// The argument of the Fock-type Airy fuctions.
      Tx x_arg;

      /// The value of the Fock-type Airy function w1.
      Tp w1_value;

      /// The derivative of the Fock-type Airy function w1.
      Tp w1_deriv;

      /// The value of the Fock-type Airy function w2.
      Tp w2_value;

      /// The derivative of the Fock-type Airy function w2.
      Tp w2_deriv;

      /// Return the Wronskian of this Fock-type Airy function state.
      constexpr Tp
      Wronskian() const noexcept
      { return w1_value * w2_deriv - w2_value * w1_deriv; }
    };

  /**
   * This struct captures the state of the cylindrical Bessel functions
   * at a given order and argument.
   */
  template<typename Tnu, typename Tx, typename Tp>
    struct cyl_bessel_t
    {
      /// The real order of the cylindrical Bessel functions.
      Tnu nu_arg;

      /// The argument of the cylindrical Bessel functions.
      Tx x_arg;

      /// The value of the Bessel function of the first kind.
      Tp J_value;

      /// The derivative of the Bessel function of the first kind.
      Tp J_deriv;

      /// The value of the Bessel function of the second kind.
      Tp N_value;

      /// The derivative of the Bessel function of the second kind.
      Tp N_deriv;

      /// Return the Wronskian of this cylindrical Bessel function state.
      constexpr Tp
      Wronskian() const noexcept
      { return J_value * N_deriv - N_value * J_deriv; }
    };

  /**
   * This struct captures the state of the Coulomb functions
   * at a given order and argument.
   */
  template<typename Teta, typename Trho, typename Tp>
    struct coulomb_t
    {
      /// The nonnegative order of the Coulomb functions.
      unsigned int l;

      /// The real parameter of the Coulomb functions.
      Teta eta_arg;

      /// The argument of the Coulomb functions.
      Trho rho_arg;

      /// The value of the regular Coulomb function.
      Tp F_value;

      /// The derivative of the regular Coulomb function.
      Tp F_deriv;

      /// The value of the irregular Coulomb function.
      Tp G_value;

      /// The derivative of the irregular Coulomb function.
      Tp G_deriv;

      /// Return the Wronskian of this Coulomb function state.
      constexpr Tp
      Wronskian() const noexcept
      { return F_value * G_deriv - G_value * F_deriv; }
    };

  /**
   * This struct captures the state of the modified cylindrical Bessel functions
   * at a given order and argument.
   */
  template<typename Tnu, typename Tx, typename Tp>
    struct cyl_mod_bessel_t
    {
      /// The real order of the modified cylindrical Bessel functions.
      Tnu nu_arg;

      /// The argument of the modified cylindrical Bessel functions.
      Tx x_arg;

      /// The value of the modified cylindrical Bessel function
      /// of the first kind.
      Tp I_value;

      /// The derivative of the modified cylindrical Bessel function
      /// of the first kind.
      Tp I_deriv;

      /// The value of the modified cylindrical Bessel function
      /// of the second kind.
      Tp K_value;

      /// The derivative of the modified cylindrical Bessel function
      /// of the second kind.
      Tp K_deriv;

      /// Return the Wronskian of this modified cylindrical Bessel function
      /// state.
      constexpr Tp
      Wronskian() const noexcept
      { return I_value * K_deriv - K_value * I_deriv; }
    };

  /**
   * Tp pretty much has to be complex.
   */
  template<typename Tnu, typename Tx, typename Tp>
    struct cyl_hankel_t
    {
      /// The real order of the cylindrical Hankel functions.
      Tnu nu_arg;

      /// The argument of the modified Hankel functions.
      Tx x_arg;

      /// The value of the cylindrical Hankel function of the first kind.
      Tp H1_value;

      /// The derivative of the cylindrical Hankel function of the first kind.
      Tp H1_deriv;

      /// The value of the cylindrical Hankel function of the second kind.
      Tp H2_value;

      /// The derivative of the cylindrical Hankel function of the second kind.
      Tp H2_deriv;

      /// Return the Wronskian of this cylindrical Hankel function state.
      constexpr Tp
      Wronskian() const noexcept
      { return H1_value * H2_deriv - H2_value * H1_deriv; }
    };

  template<typename Tn, typename Tx, typename Tp>
    struct sph_bessel_t
    {
      /// The integral order of the spherical Bessel functions.
      Tn n_arg;

      /// The argument of the spherical Bessel functions.
      Tx x_arg;

      /// The value of the spherical Bessel function of the first kind.
      Tp j_value;

      /// The derivative of the spherical Bessel function of the first kind.
      Tp j_deriv;

      /// The value of the spherical Bessel function of the second kind.
      Tp n_value;

      /// The derivative of the spherical Bessel function of the second kind.
      Tp n_deriv;

      /// Return the Wronskian of this spherical Bessel function state.
      constexpr Tp
      Wronskian() const noexcept
      { return j_value * n_deriv - n_value * j_deriv; }
    };

  template<typename Tn, typename Tx, typename Tp>
    struct sph_mod_bessel_t
    {
      /// The integral order of the modified spherical Bessel functions.
      Tn n_arg;

      /// The argument of the modified spherical Bessel functions.
      Tx x_arg;

      /// The value of the modified spherical Bessel function
      /// of the first kind.
      Tp i_value;

      /// The derivative of the modified spherical Bessel function
      /// of the first kind.
      Tp i_deriv;

      /// The value of the modified spherical Bessel function
      /// of the second kind.
      Tp k_value;

      /// The derivative of the modified spherical Bessel function
      /// of the second kind.
      Tp k_deriv;

      /// Return the Wronskian of this modified cylindrical Bessel function
      /// state.
      constexpr Tp
      Wronskian() const noexcept
      { return i_value * k_deriv - k_value * i_deriv; }
    };

  /**
   * Tp pretty much has to be complex.
   */
  template<typename Tn, typename Tx, typename Tp>
    struct sph_hankel_t
    {
      /// The integral order of the spherical Hankel functions.
      Tn n_arg;

      /// The argument of the spherical Hankel functions.
      Tx x_arg;

      /// The velue of the spherical Hankel function of the first kind.
      Tp h1_value;

      /// The derivative of the spherical Hankel function of the first kind.
      Tp h1_deriv;

      /// The velue of the spherical Hankel function of the second kind.
      Tp h2_value;

      /// The derivative of the spherical Hankel function of the second kind.
      Tp h2_deriv;

      /// Return the Wronskian of this cylindrical Hankel function state.
      constexpr Tp
      Wronskian() const noexcept
      { return h1_value * h2_deriv - h2_value * h1_deriv; }
    };

  template<typename Tp>
    struct gappa_pq_t
    {
      /// 
      Tp gappa_p_value;

      /// 
      Tp gappa_q_value;
    };

  /**
   * The log of the absolute value of the gamma function
   * The sign of the exponentiated log(gamma) is stored in sign.
   */
  template<typename Tp>
    struct lgamma_t
    {
      /// The value log gamma function.
      Tp lgamma_value;

      /// The sign of the exponent of the log gamma value.
      int lgamma_sign;
    };

  /**
   * The sign of the exponentiated log(gamma) is appied to the tgamma value.
   */
  template<typename Tp>
    struct gamma_inc_t
    {
      /// The value of the total gamma function.
      Tp tgamma_value;
      /// The value of the log of the incomplete gamma function
      Tp lgamma_value;
    };

  /**
   * @brief A structure for the gamma functions required by the Temme series
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
   * The accuracy requirements on this are high for @f$ |\mu| < 0 @f$.
   */
  template<typename Tp>
    struct gamma_temme_t
    {
      /// The input parameter of the gamma functions
      Tp mu_arg;

      /// The output function @f$ 1/\Gamma(1 + \mu) @f$
      Tp gamma_plus_value;

      /// The output function @f$ 1/\Gamma(1 - \mu) @f$
      Tp gamma_minus_value;

      /// The output function @f$ \Gamma_1(\mu) @f$
      Tp gamma_1_value;

      /// The output function @f$ \Gamma_2(\mu) @f$
      Tp gamma_2_value;
    };

  /**
   * @brief A structure for Stirling numbers of the first kind.
   */
  template<typename Tp>
    struct stirling_1_t
    {
      using iterator = typename std::vector<Tp>::iterator;
      using const_iterator = typename std::vector<Tp>::const_iterator;

      std::vector<Tp> sigma;

      unsigned int
      degree() const noexcept
      { return sigma.size() - 1; }

      Tp
      operator[](unsigned int k) const noexcept
      { return k < sigma.size() ? sigma[k] : Tp{0}; }

      template<typename Up>
	auto
	operator()(Up x) const noexcept
	{
	  const auto n = sigma.size() - 1;
	  auto poly = sigma[n];
	  for (unsigned int i = 1; i < n; ++i)
	    poly = sigma[n - i] + x * poly;
	}

      iterator
      begin() noexcept
      { return sigma.begin(); }

      iterator
      end() noexcept
      { return sigma.end(); }

      const_iterator
      begin() const noexcept
      { return this->sigma.begin(); }

      const_iterator
      end() const noexcept
      { return this->sigma.end(); }
    };

  /**
   * @brief A structure for Stirling numbers of the first kind.
   */
  template<typename Tp>
    struct stirling_2_t
    {
      using iterator = typename std::vector<Tp>::iterator;
      using const_iterator = typename std::vector<Tp>::const_iterator;

      std::vector<Tp> S;

      unsigned int
      degree() const noexcept
      { return S.size() - 1; }

      Tp
      operator[](unsigned int k) const noexcept
      { return k < S.size() ? S[k] : Tp{0}; }

      /// Return the Bell polynomial.
      template<typename Up>
	auto
	operator()(Up x) const noexcept
	{
	  const auto n = S.size() - 1;
	  auto poly = S[n];
	  for (unsigned int i = 1; i < n; ++i)
	    poly = S[n - i] + x * poly;
	}

      iterator
      begin() noexcept
      { return S.begin(); }

      iterator
      end() noexcept
      { return S.end(); }

      const_iterator
      begin() const noexcept
      { return this->S.begin(); }

      const_iterator
      end() const noexcept
      { return this->S.end(); }
    };

} // namespace emsr

#endif // SPECFUN_STATE_H
