#ifndef SF_DISTRIBUTIONS_H
#define SF_DISTRIBUTIONS_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_distributions.tcc>

namespace emsr
{

  /**
   * @brief Return the gamma propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    gamma_pdf(Ta alpha, Tb beta, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::gamma_pdf<type>(alpha, beta, x);
    }

  /**
   * @brief Return the gamma cumulative propability distribution function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    gamma_p(Ta alpha, Tb beta, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::gamma_p<type>(alpha, beta, x);
    }
   */

  /**
   * @brief Return the normal probability density function.
   *
   * The formula for the normal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(x-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename Tmu, typename Tsig, typename Tp>
    inline emsr::fp_promote_t<Tmu, Tsig, Tp>
    normal_pdf(Tmu mu, Tsig sigma, Tp x)
    {
      using type = emsr::fp_promote_t<Tmu, Tsig, Tp>;
      return emsr::detail::normal_pdf<type>(mu, sigma, x);
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
  template<typename Tmu, typename Tsig, typename Tp>
    inline emsr::fp_promote_t<Tmu, Tsig, Tp>
    normal_p(Tmu mu, Tsig sigma, Tp x)
    {
      using type = emsr::fp_promote_t<Tmu, Tsig, Tp>;
      return emsr::detail::normal_p<type>(mu, sigma, x);
    }

  /**
   * @brief Return the lognormal probability density function.
   *
   * The formula for the lognormal probability density function is
   * @f[
   *   f(x|\mu,\sigma) = \frac{e^{(\ln{x}-\mu)^2/2\sigma^2}}{\sigma\sqrt{2\pi}}
   * @f]
   */
  template<typename Tmu, typename Tsig, typename Tp>
    inline emsr::fp_promote_t<Tmu, Tsig, Tp>
    lognormal_pdf(Tmu mu, Tsig sigma, Tp x)
    {
      using type = emsr::fp_promote_t<Tmu, Tsig, Tp>;
      return emsr::detail::lognormal_pdf<type>(mu, sigma, x);
    }

  /**
   * @brief Return the lognormal cumulative probability density function.
   *
   * The formula for the lognormal cumulative probability density function is
   * @f[
   *   F(x|\mu,\sigma)
   *     = \frac{1}{2}\left[1-erf\left(
   *                  \frac{\ln{x}-\mu}{\sqrt{2}\sigma}\right)\right]
   * @f]
   */
  template<typename Tmu, typename Tsig, typename Tp>
    inline emsr::fp_promote_t<Tmu, Tsig, Tp>
    lognormal_p(Tmu mu, Tsig sigma, Tp x)
    {
      using type = emsr::fp_promote_t<Tmu, Tsig, Tp>;
      return emsr::detail::lognormal_p<type>(mu, sigma, x);
    }

  /**
   * @brief Return the exponential probability density function.
   *
   * The formula for the exponential probability density function is
   * @f[
   *   f(x|\lambda) = \lambda e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tlam, typename Tp>
    inline emsr::fp_promote_t<_Tlam, Tp>
    exponential_pdf(_Tlam lambda, Tp x)
    {
      using type = emsr::fp_promote_t<_Tlam, Tp>;
      return emsr::detail::exponential_pdf<type>(lambda, x);
    }

  /**
   * @brief Return the exponential cumulative probability density function.
   *
   * The formula for the exponential cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-\lambda x} \mbox{ for } x >= 0
   * @f]
   */
  template<typename _Tlam, typename Tp>
    inline emsr::fp_promote_t<_Tlam, Tp>
    exponential_p(_Tlam lambda, Tp x)
    {
      using type = emsr::fp_promote_t<_Tlam, Tp>;
      return emsr::detail::exponential_p<type>(lambda, x);
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
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    weibull_pdf(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::weibull_pdf<type>(a, b, x);
    }

  /**
   * @brief Return the Weibull cumulative probability density function.
   *
   * The formula for the Weibull cumulative probability density function is
   * @f[
   *   F(x|\lambda) = 1 - e^{-(x / b)^a} \mbox{ for } x >= 0
   * @f]
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    weibull_p(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::weibull_p<type>(a, b, x);
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
  template<typename _Tt, typename Tp>
    emsr::fp_promote_t<Tp>
    student_t_pdf(_Tt t, unsigned int nu)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::student_t_pdf<type>(t, nu);
    }

  /**
   * @brief  Return the Students T probability function.
   *
   * The students T propability function is related to the incomplete beta function:
   * @f[
   *   A(t|\nu) = 1 - I_{\frac{\nu}{\nu + t^2}}(\frac{\nu}{2}, \frac{1}{2})
   * @f]
   * (see ibeta).
   *
   * @param t 
   * @param nu 
   */
  template<typename _Tt, typename Tp>
    emsr::fp_promote_t<Tp>
    student_t_p(_Tt t, unsigned int nu)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::student_t_p<type>(t, nu);
    }

  /**
   * @brief  Return the F-distribution propability function.
   * This returns the probability that the observed chi-square
   * for a correct model exceeds the value @f$ \chi^2 @f$.
   *
   * The f-distribution propability function is related
   * to the incomplete beta function:
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
    emsr::fp_promote_t<Tp>
    fisher_f_pdf(Tp F, unsigned int nu1, unsigned int nu2)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::fisher_f_pdf<type>(F, nu1, nu2);
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
    emsr::fp_promote_t<Tp>
    fisher_f_p(Tp F, unsigned int nu1, unsigned int nu2)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::fisher_f_p<type>(F, nu1, nu2);
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
    emsr::fp_promote_t<Tp>
    binomial_pdf(Tp p, unsigned int n, unsigned int k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::binomial_pdf<type>(p, n, k);
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
    emsr::fp_promote_t<Tp>
    binomial_p(Tp p, unsigned int n, unsigned int k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::binomial_p<type>(p, n, k);
    }

  /**
   * @brief  Return the logistic probability density function.
   *
   * The formula for the logistic probability density function is
   * @f[
   *     f(x| a, b) = \frac{e^{(x - a)/b}}{b[1 + e^{(x - a)/b}]^2}
   * @f]
   * where @f$b > 0@f$.
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    logistic_pdf(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::logistic_pdf<type>(a, b, x);
    }

  /**
   * @brief  Return the logistic cumulative distribution function.
   *
   * The formula for the logistic probability function is
   * @f[
   *     P(x| a, b) = \frac{e^{(x - a)/b}}{1 + e^{(x - a)/b}}
   * @f]
   * where @f$b > 0@f$.
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    logistic_p(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::logistic_p<type>(a, b, x);
    }

} // namespace emsr

#endif // SF_DISTRIBUTIONS_H
