
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

/** @file bits/sf_distributions.tcc
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
//     Section 6, pp. 253-266
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Numerical Recipes in C, by W. H. Press, S. A. Teukolsky,
//     W. T. Vetterling, B. P. Flannery, Cambridge University Press (1992),
//     2nd ed, pp. 213-216
// (4) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.

#ifndef SF_DISTRIBUTIONS_TCC
#define SF_DISTRIBUTIONS_TCC 1

#include <stdexcept>

namespace emsr
{
namespace detail
{

  /**
   *  @brief  Return the chi-squared propability function.
   *  This returns the probability that the observed chi-squared for a correct model
   *  is less than the value @f$ \chi^2 @f$.
   *
   *  The chi-squared propability function is related
   *  to the normalized lower incomplete gamma function:
   *  @f[
   *    P(\chi^2|\nu) = \Gamma_P(\frac{\nu}{2}, \frac{\chi^2}{2})
   *  @f]
   */
  template<typename Tp>
    Tp
    chi_squared_pdf(Tp chi2, unsigned int nu)
    {
      if (std::isnan(chi2))
	return emsr::quiet_NaN(chi2);
      else if (chi2 < Tp{0})
	throw std::domain_error("chi_squared_p: chi-squared is negative");
      else
	return gamma_p(Tp(nu) / Tp{2}, chi2 / Tp{2});
    }

  /**
   *  @brief  Return the complementary chi-squared propability function.
   *  This returns the probability that the observed chi-squared for a correct model
   *  is greater than the value @f$ \chi^2 @f$.
   *
   *  The complementary chi-squared propability function is related
   *  to the normalized upper incomplete gamma function:
   *  @f[
   *    Q(\chi^2|\nu) = \Gamma_Q(\frac{\nu}{2}, \frac{\chi^2}{2})
   *  @f]
   */
  template<typename Tp>
    Tp
    chi_squared_pdfc(Tp chi2, unsigned int nu)
    {
      if (std::isnan(chi2) || std::isnan(nu))
	return emsr::quiet_NaN(chi2);
      else if (chi2 < Tp{0})
	throw std::domain_error("chi_square_pdfc: chi-squared is negative");
      else
	return gamma_q(Tp(nu) / Tp{2}, chi2 / Tp{2});
    }

  /**
   * @brief Return the gamma propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename Tp>
    Tp
    gamma_pdf(Tp alpha, Tp beta, Tp x)
    {
      if (std::isnan(alpha) || std::isnan(beta) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return std::pow(beta, alpha) * std::pow(x, alpha - Tp{1})
	     * std::exp(beta * x) / gamma(alpha);
    }

  /**
   * @brief Return the gamma cumulative propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename Tp>
    Tp
    gamma_p(Tp alpha, Tp beta, Tp x)
    {
      if (std::isnan(alpha) || std::isnan(beta) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return tgamma_lower(alpha, beta * x)
	     / gamma(alpha);
    }

  /**
   * @brief Return the gamma complementary cumulative propability
   *        distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                            (x/\beta)^{\alpha - 1} e^{-x/\beta} 
   * @f]
   */
  template<typename Tp>
    Tp
    gamma_q(Tp alpha, Tp beta, Tp x)
    {
      if (std::isnan(alpha) || std::isnan(beta) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return tgamma(alpha, beta * x)
	     / gamma(alpha);
    }


  /**
   * @brief Return the Rice probability density function.
   *
   * The formula for the Rice probability density function is
   * @f[
   *   p(x|\nu,\sigma) = \frac{x}{\sigma^2}
   *                     \exp\left(-\frac{x^2+\nu^2}{2\sigma^2}\right)
   *                     I_0\left(\frac{x \nu}{\sigma^2}\right)
   * @f]
   * where @f$I_0(x)@f$ is the modified Bessel function of the first kind
   * of order 0 and @f$\nu >= 0@f$ and @f$\sigma > 0@f$.
   */
  template<typename Tp>
    Tp
    rice_pdf(Tp nu, Tp sigma, Tp x)
    {
      if (std::isnan(nu) || std::isnan(sigma))
	return emsr::quiet_NaN(x);
      else
	{
	  auto sigma2 = sigma * sigma;
	  return (x / sigma2)
               * std::exp(-(x * x + nu * nu) / (Tp{2} * sigma2))
               * cyl_bessel_i(Tp{0}, (x * nu) / (sigma2));
	}
    }


  /**
   * @brief Return the normal probability density function.
   *
   * The formula for the normal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(x-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename Tp>
    Tp
    normal_pdf(Tp mu, Tp sigma, Tp x)
    {
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
      const auto s_sqrt_2pi = s_sqrt_2 * s_sqrt_pi;
      if (std::isnan(mu) || std::isnan(sigma))
	return emsr::quiet_NaN(x);
      else
	{
	  x -= mu;
	  x /= sigma;
	  x *= x;
	  x /= Tp{2};
	  return std::exp(-x) / (sigma * s_sqrt_2pi);
	}
    }

  /**
   * @brief Return the normal cumulative probability density function.
   *
   * The formula for the normal cumulative probability density function is
   * @f[
   *     F(x|\mu,\sigma)
   *        = \frac{1}{2}\left[ 1-erf(\frac{x-\mu}{\sqrt{2}\sigma}) \right]
   * @f]
   */
  template<typename Tp>
    Tp
    normal_p(Tp mu, Tp sigma, Tp x)
    {
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      if (std::isnan(mu) || std::isnan(sigma) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else
	return Tp{0.5L} * (Tp{1} + std::erf((x - mu) / (sigma * s_sqrt_2)));
    }


  /**
   * @brief Return the lognormal probability density function.
   *
   * The formula for the lognormal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(\ln{x}-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename Tp>
    Tp
    lognormal_pdf(Tp nu, Tp sigma, Tp x)
    {
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      const auto s_sqrt_pi = emsr::sqrtpi_v<Tp>;
      const auto s_sqrt_2pi = s_sqrt_2 * s_sqrt_pi;
      if (std::isnan(nu) || std::isnan(sigma) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
        throw std::domain_error("lognormal_pdf: argument x must be positive.");
      else if (x == Tp{0})
        return Tp{0};
      else
	{
          auto arg = std::log(x);
	  arg -= nu;
	  arg /= sigma;
	  arg *= arg;
	  arg /= Tp{2};
	  return std::exp(-arg) / (sigma * s_sqrt_2pi);
	}
    }

  /**
   * @brief Return the lognormal cumulative probability density function.
   *
   * The formula for the lognormal cumulative probability density function is
   * @f[
   *   F(x|\mu,\sigma)
   *     = \frac{1}{2}\left[ 1-erf(\frac{\ln{x}-\mu}{\sqrt{2}\sigma}) \right]
   * @f]
   */
  template<typename Tp>
    Tp
    lognormal_p(Tp mu, Tp sigma, Tp x)
    {
      const auto s_sqrt_2 = emsr::sqrt2_v<Tp>;
      if (std::isnan(mu) || std::isnan(sigma) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
        throw std::domain_error("lognormal_p: argument x must be positive.");
      else if (x == Tp{0})
        return Tp{0};
      else
	return Tp{0.5L} * (Tp{1} + std::erf((std::log(x) - mu)
					    / (sigma * s_sqrt_2)));
    }


  /**
   * @brief Return the exponential probability density function.
   *
   * The formula for the exponential probability density function is
   * @f[
   *   f(x|\lambda) = \lambda e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename Tp>
    Tp
    exponential_pdf(Tp lambda, Tp x)
    {
      if (std::isnan(lambda) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return Tp{0};
      else
	return lambda * std::exp(-lambda * x);
    }

  /**
   * @brief Return the exponential cumulative probability density function.
   *
   * The formula for the exponential cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename Tp>
    Tp
    exponential_p(Tp lambda, Tp x)
    {
      if (std::isnan(lambda) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return Tp{0};
      else
	return Tp{1} - std::exp(-lambda * x);
    }

  /**
   * @brief Return the complement of the exponential cumulative
   * probability density function.
   *
   * The formula for the complement of the exponential cumulative
   * probability density function is
   * @f[
   *   F(x|\lambda) = e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename Tp>
    Tp
    exponential_q(Tp lambda, Tp x)
    {
      if (std::isnan(lambda) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return Tp{0};
      else
	return std::exp(-lambda * x);
    }


  /**
   * @brief Return the Weibull probability density function.
   *
   * The formula for the Weibull probability density function is
   * @f[
   *   f(x | a, b) = \frac{a}{b}
   *                 \left(\frac{x}{b} \right)^{a-1}
   *                 \exp{-\left(\frac{x}{b}\right)^a}
   *                 \mbox{ for } x >= 0
   * @f]
   */
  template<typename Tp>
    Tp
    weibull_pdf(Tp a, Tp b, Tp x)
    {
      if (std::isnan(a) || std::isnan(b) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return Tp{0};
      else
	return (a / b) * std::pow(x / b, a - Tp{1})
 			   * std::exp(-std::pow(x / b, a));
    }

  /**
   * @brief Return the Weibull cumulative probability density function.
   *
   * The formula for the Weibull cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-(x / b)^a} \mbox{ for } x >= 0
   * @f]
   */
  template<typename Tp>
    Tp
    weibull_p(Tp a, Tp b, Tp x)
    {
      if (std::isnan(a) || std::isnan(b) || std::isnan(x))
	return emsr::quiet_NaN(x);
      else if (x < Tp{0})
	return Tp{0};
      else
	return Tp{1} - std::exp(-std::pow(x / b, a));
    }

  /**
   * @brief  Return the Students T probability density.
   *
   * The students T propability density is:
   * @f[
   *   A(t|\nu) = 1 - I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   *   A(t|\nu) = 
   * @f]
   *
   * @param t 
   * @param nu 
   */
  template<typename Tp>
    Tp
    student_t_pdf(Tp t, unsigned int nu)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      if (std::isnan(t))
	return emsr::quiet_NaN(t);
      else
	return gamma(Tp(nu + 1) / Tp{2})
	     * std::pow((Tp(nu) + t * t) / nu, -Tp(nu + 1) / 2 )
	     / gamma(Tp(nu) / Tp{2}) / std::sqrt(Tp(nu) * s_pi);
    }

  /**
   * @brief  Return the Students T probability function.
   *
   * The students T propability function is related to the incomplete beta function:
   * @f[
   *   A(t|\nu) = 1 - I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   *   A(t|\nu) = 
   * @f]
   *
   * @param t 
   * @param nu 
   */
  template<typename Tp>
    Tp
    student_t_p(Tp t, unsigned int nu)
    {
      if (std::isnan(t))
	return emsr::quiet_NaN(t);
      else
	return beta_inc(Tp{0.5L}, Tp(nu) / Tp{2},
			  t * t / (Tp(nu) + t * t));
    }

  /**
   * @brief  Return the complement of the Students T probability function.
   *
   * The complement of the students T propability function is:
   * @f[
   *   A_c(t|\nu) = I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   * 		  = 1 - A(t|\nu)
   * @f]
   *
   * @param t 
   * @param nu 
   */
  template<typename Tp>
    Tp
    student_t_q(Tp t, unsigned int nu)
    {
      if (std::isnan(t))
	return emsr::quiet_NaN(t);
      else
	return beta_inc(Tp(nu) / Tp{2}, Tp{0.5L},
			  Tp(nu) / (Tp(nu) + t * t));
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square for a correct model
   * exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related to the incomplete beta function:
   * @f[
   *   Q(F|\nu_1, \nu_2) = I_{\frac{\nu_2}{\nu_2 + \nu_1 F}}
   * 			     (\frac{\nu_2}{2}, \frac{\nu_1}{2})
   * @f]
   *
   * @param nu1 The number of degrees of freedom of sample 1
   * @param nu2 The number of degrees of freedom of sample 2
   * @param F The F statistic
   */
  template<typename Tp>
    Tp
    fisher_f_pdf(Tp F, unsigned int nu1, unsigned int nu2)
    {
      if (std::isnan(F))
	return emsr::quiet_NaN(F);
      else if (F < Tp{0})
	throw std::domain_error("f_p: F is negative");
      else
	return std::sqrt(std::pow(Tp(nu1) * F, Tp(nu1))
		       * std::pow(Tp(nu2), Tp(nu2))
		/ std::pow(Tp(nu1) * F + Tp(nu2), Tp(nu1 + nu2)))
	    / F / beta(Tp(nu1) / Tp{2}, Tp(nu2) / Tp{2});
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square for a correct model
   * exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related to the incomplete beta function:
   * @f[
   *   Q(F|\nu_1, \nu_2) = I_{\frac{\nu_2}{\nu_2 + \nu_1 F}}
   * 			     (\frac{\nu_2}{2}, \frac{\nu_1}{2})
   * @f]
   *
   * @param nu1 The number of degrees of freedom of sample 1
   * @param nu2 The number of degrees of freedom of sample 2
   * @param F The F statistic
   */
  template<typename Tp>
    Tp
    fisher_f_p(Tp F, unsigned int nu1, unsigned int nu2)
    {
      if (std::isnan(F))
	return emsr::quiet_NaN(F);
      else if (F < Tp{0})
	throw std::domain_error("f_p: F is negative");
      else
	return beta_inc(Tp(nu2) / Tp{2}, Tp(nu1) / Tp{2},
			  Tp(nu2) / (Tp(nu2) + nu1 * F));
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square for a correct model
   * exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related to the incomplete beta function:
   * @f[
   *   P(F|\nu_1, \nu_2) = 1 - I_{\frac{\nu_2}{\nu_2 + \nu_1 F}}
   * 			     (\frac{\nu_2}{2}, \frac{\nu_1}{2})
   * 			 = 1 - Q(F|\nu_1, \nu_2)
   * @f]
   *
   * @param F 
   * @param nu1 
   * @param nu2 
   */
  template<typename Tp>
    Tp
    fisher_f_q(Tp F, unsigned int nu1, unsigned int nu2)
    {
      if (std::isnan(F))
	return emsr::quiet_NaN(F);
      else if (F < Tp{0})
	throw std::domain_error("f_q: F is negative");
      else
	return beta_inc(Tp(nu1) / Tp{2}, Tp(nu2) / Tp{2},
			  nu1 * F / (Tp(nu2) + nu1 * F));
    }

  /**
   * @brief  Return the binomial probability mass function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   f(k|n,p) = \binom{n}{k}p^k(1-p)^{n-k}
   * @f]
   *
   * @param p 
   * @param n 
   * @param k 
   */
  template<typename Tp>
    Tp
    binomial_pdf(Tp p, unsigned int n, unsigned int k)
    {
      if (std::isnan(p))
	return emsr::quiet_NaN(p);
      else if (p < Tp{0} || p > Tp{1})
	throw std::domain_error("binomial_p: probability is out of range");
      else if (k > n)
	return Tp{0};
      else if (n == 0)
	return Tp{1};
      else if (k == 0)
	return std::pow(Tp{1} - p, n);
      else if (k == n)
	return std::pow(p, n);
      else
	return binomial<Tp>(n, k)
	     * std::pow(p, k)
	     * std::pow(Tp{1} - p, n - k);
    }

  /**
   * @brief  Return the binomial cumulative distribution function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   P(k|n,p) = I_p(k, n-k+1)
   * @f]
   *
   * @param p 
   * @param n 
   * @param k 
   */
  template<typename Tp>
    Tp
    binomial_p(Tp p, unsigned int n, unsigned int k)
    {
      if (std::isnan(p))
	return emsr::quiet_NaN(p);
      else if (p < Tp{0} || p > Tp{1})
	throw std::domain_error("binomial_p: probability is out of range");
      else if (k == 0)
	return Tp{1};
      else if (k > n)
	return Tp{0};
      else
	return beta_inc(Tp(k), Tp(n - k - 1), p);
    }

  /**
   * @brief  Return the complementary binomial cumulative distribution function.
   *
   * The binomial cumulative distribution function is related
   * to the incomplete beta function:
   * @f[
   *   Q(k|n,p) = I_{1-p}(n-k+1, k)
   * @f]
   *
   * @param p 
   * @param n 
   * @param k 
   */
  template<typename Tp>
    Tp
    binomial_q(Tp p, unsigned int n, unsigned int k)
    {
      if (std::isnan(p))
	return emsr::quiet_NaN(p);
      else if (p < Tp{0} || p > Tp{1})
	throw std::domain_error("binomial_q: probability is out of range");
      else if (k == 0)
	return Tp{1};
      else if (k > n)
	return Tp{0};
      else
	return beta_inc(Tp(n - k - 1), Tp(k), Tp{1} - p);
    }

  /**
   * @brief  Return the logistic probability density function.
   *
   * The formula for the logistic probability density function is
   * @f[
   *     p(x| mu, s) = \frac{e^{-(x - \mu)/s}}{s[1 + e^{-(x - \mu)/s}]^2}
   * @f]
   * where @f$s > 0@f$.
   */
  template<typename Tp>
    Tp
    logistic_pdf(Tp mu, Tp s, Tp x)
    {
      const auto arg = -(x - mu) / s;
      const auto exparg = std::exp(arg);
      return exparg / (s * (Tp{1} + exparg) * (Tp{1} + exparg));
    }

  /**
   * @brief  Return the logistic cumulative distribution function.
   *
   * The formula for the logistic probability function is
   * @f[
   *     cdf(x| \mu, s) = \frac{1}{1 + e^{-(x - \mu)/s}}
   * @f]
   * where @f$s > 0@f$.
   */
  template<typename Tp>
    Tp
    logistic_p(Tp mu, Tp s, Tp x)
    {
      const auto arg = -(x - mu) / s;
      const auto exparg = std::exp(arg);
      return exparg / (Tp{1} + exparg);
    }

  template<typename Tp>
    Tp
    cauchy_p(Tp a, Tp b, Tp x)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      return Tp{0.5L} + std::atan((x - a) / b) / s_pi;
    }

  template<typename Tp>
    Tp
    beta_p(Tp a, Tp b, Tp x)
    {
      if (x < Tp{0})
	return Tp{0};
      else if (x > Tp{1})
	return Tp{1};
      else
	return beta_inc(a, b, x);
    }

  /**
   * @f[
   *     P(K <= x) = 1 - e^{-2x^2} + e^{-2 \cdot 4 x^2}
   *            + e^{-2 \cdot 9 x^2} - e^{-2 \cdot 16 x^2} + ...
   * @f]
   */
  template<typename Tp>
    Tp
    kolmogorov_p(Tp a, Tp b, Tp x)
    {
      return Tp{1} - std::exp(-Tp{2} * x * x);
    }

} // namespace detail
} // namespace emsr

#endif // SF_DISTRIBUTIONS_TCC
