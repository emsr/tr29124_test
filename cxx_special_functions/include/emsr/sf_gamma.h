#ifndef SF_GAMMA_H
#define SF_GAMMA_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_gamma.tcc>

namespace emsr
{

  // Factorial

  /**
   * @brief Return the factorial @f$ n! @f$ of the argument as a real number.
   * @f[
   *   n! = 1 \times 2 \times ... \times n, 0! = 1
   * @f]
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    factorial(unsigned int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::factorial<type>(n);
    }

  // Double factorial

  /**
   * @brief Return the double factorial @f$ n!! @f$ of the argument
   * as a real number.
   * @f[
   *   n!! = n(n-2)...(2), 0!! = 1
   * @f]
   * for even @c n and
   * @f[
   *   n!! = n(n-2)...(1), (-1)!! = 1
   * @f]
   * for odd @c n.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    double_factorial(int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::double_factorial<type>(n);
    }

  // Log factorial

  /**
   * @brief Return the logarithm of the factorial @f$ ln(n!) @f$ of the argument
   * as a real number.
   * @f[
   *   n! = 1 \times 2 \times ... \times n, 0! = 1
   * @f]
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    lfactorial(unsigned int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::log_factorial<type>(n);
    }

  // Log double factorial

  /**
   * @brief Return the logarithm of the double factorial @f$ ln(n!!) @f$
   * of the argument as a real number.
   * @f[
   *   n!! = n(n-2)...(2), 0!! = 1
   * @f]
   * for even @c n and
   * @f[
   *   n!! = n(n-2)...(1), (-1)!! = 1
   * @f]
   * for odd @c n.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    ldouble_factorial(int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::log_double_factorial<type>(n);
    }

  //  Log rising factorials

  /**
   * @brief  Return the logarithm of the rising factorial function
   * or the (upper) Pochhammer symbol.
   * The rising factorial function is defined for integer order by
   * @f[
   *   a^{\overline{\nu}} = \Gamma(a + \nu) / \Gamma(n)
   *	     = \prod_{k=0}^{\nu-1} (a + k), a^{\overline{0}} \equiv 1
   * @f]
   * Thus this function returns
   * @f[
   *   ln[a^{\overline{\nu}}] = ln[\Gamma(a + \nu)] - ln[\Gamma(\nu)],
   *      ln[a^{\overline{0}}] \equiv 0
   * @f]
   * Many notations exist for this function: @f$ (a)_\nu @f$, called
   * the Pochhammer function (esp. in the literature of special functions), and
   * @f[ \left[ \begin{array}{c}
   *	  a \\
   *	  \nu \end{array} \right]
   * @f], and others.
   */
  template<typename Tp, typename Tnu>
    inline emsr::fp_promote_t<Tp, Tnu>
    lrising_factorial(Tp a, Tnu nu)
    {
      using type = emsr::fp_promote_t<Tp, Tnu>;
      return emsr::detail::log_rising_factorial<type>(a, nu);
    }

  //  Log falling factorials

  /**
   * @brief  Return the logarithm of the falling factorial function
   * or the lower Pochhammer symbol.
   * The falling factorial function is defined by
   * @f[
   *   a^{\underline{n}} = \frac{\Gamma(a + 1)}{\Gamma(a - \nu + 1)}
   *	     = \prod_{k=0}^{n-1} (a - k)
   * @f]
   * where @f$ a^{\underline{0}} \equiv 1 @f$.
   * In particular, @f$ n^{\underline{n}} = n! @f$.
   * Thus this function returns
   * @f[
   *   ln[a^{\underline{n}}] = ln[\Gamma(a + 1)] - ln[\Gamma(a - \nu + 1)]
   * @f]
   * where @f$ ln[a^{\underline{0}}] \equiv 0 @f$.
   * Many notations exist for this function: @f$ (a)_\nu @f$,
   * @f[
   *    \{ \begin{array}{c}
   *	  a \\
   *	  \nu \end{array} \}
   * @f], and others.
   */
  template<typename Tp, typename Tnu>
    inline emsr::fp_promote_t<Tp, Tnu>
    lfalling_factorial(Tp a, Tnu nu)
    {
      using type = emsr::fp_promote_t<Tp, Tnu>;
      return emsr::detail::log_falling_factorial<type>(a, nu);
    }

  //  Rising factorials

  /**
   * @brief  Return the rising factorial function
   * or the (upper) Pochhammer function.
   * The rising factorial function is defined by
   * @f[
   *   a^{\overline{\nu}} = \Gamma(a + \nu) / \Gamma(\nu)
   * @f]
   * Many notations exist for this function: @f$ (a)_\nu @f$, called
   * the Pochhammer function (esp. in the literature of special functions), and
   * @f[ \left[ \begin{array}{c}
   *	  a \\
   *	  \nu \end{array} \right]
   * @f], and others.
   */
  template<typename Tp, typename Tnu>
    inline emsr::fp_promote_t<Tp, Tnu>
    rising_factorial(Tp a, Tnu nu)
    {
      using type = emsr::fp_promote_t<Tp, Tnu>;
      return emsr::detail::rising_factorial<type>(a, nu);
    }

  //  Falling factorials

  /**
   * @brief  Return the falling factorial function
   * or the lower Pochhammer symbol for real argument @c a
   * and integral order @c n.
   * The falling factorial function is defined by
   * @f[
   *   a^{\underline{n}} = \prod_{k=0}^{n-1} (a - k)
   *	     = \Gamma(a + 1) / \Gamma(a - n + 1)
   * @f]
   * where @f$ a^{\underline{0}} \equiv 1 @f$.
   * In particular, @f$ n^{\underline{n}} = n! @f$.
   */
  template<typename Tp, typename Tnu>
    inline emsr::fp_promote_t<Tp, Tnu>
    falling_factorial(Tp a, Tnu nu)
    {
      using type = emsr::fp_promote_t<Tp, Tnu>;
      return emsr::detail::falling_factorial<type>(a, nu);
    }

  // Binomial coefficient

  /**
   * @brief Return the binomial coefficient as a real number.
   * The binomial coefficient is given by:
   * @f[
   *   \binom{n}{k} = \frac{n!}{(n-k)! k!}
   * @f]
   * The binomial coefficients are generated by:
   * @f[
   *   \left(1 + t\right)^n = \sum_{k=0}^n \binom{n}{k} t^k
   * @f]
   *
   * @param n The first argument of the binomial coefficient.
   * @param k The second argument of the binomial coefficient.
   * @return  The binomial coefficient.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    binomial(unsigned int n, unsigned int k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::binomial<type>(n, k);
    }

  // Log binomial coefficient

  /**
   * @brief Return the logarithm of the binomial coefficient as a real number.
   * The binomial coefficient is given by:
   * @f[
   *   \binom{n}{k} = \frac{n!}{(n-k)! k!}
   * @f]
   * The binomial coefficients are generated by:
   * @f[
   *   \left(1 + t\right)^n = \sum_{k=0}^n \binom{n}{k} t^k
   * @f]
   *
   * @param n The first argument of the binomial coefficient.
   * @param k The second argument of the binomial coefficient.
   * @return  The logarithm of the binomial coefficient.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    lbinomial(unsigned int n, unsigned int k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::log_binomial<type>(n, k);
    }

  // Log Gamma function for real argument.

  /**
   * Return the logarithm of the gamma function for real argument.
   */
  template<typename _Ta>
    inline emsr::fp_promote_t<_Ta>
    lgamma(_Ta a)
    {
      using type = emsr::fp_promote_t<_Ta>;
      return emsr::detail::log_gamma<type>(a);
    }

  // Log Gamma function for complex argument.

  /**
   * Return the logarithm of the gamma function for complex argument.
   */
  template<typename _Ta>
    inline std::complex<emsr::fp_promote_t<_Ta>>
    lgamma(std::complex<_Ta> a)
    {
      using type = std::complex<emsr::fp_promote_t<_Ta>>;
      return emsr::detail::log_gamma<type>(a);
    }

  // Gamma function for real argument.

  /**
   * Return the gamma function for real argument.
   */
  template<typename _Ta>
    inline emsr::fp_promote_t<_Ta>
    tgamma(_Ta a)
    {
      using type = emsr::fp_promote_t<_Ta>;
      return emsr::detail::gamma<type>(a);
    }

  // Gamma function for complex argument.

  /**
   * Return the gamma function for complex argument.
   */
  template<typename _Ta>
    inline std::complex<emsr::fp_promote_t<_Ta>>
    tgamma(std::complex<_Ta> a)
    {
      using type = std::complex<emsr::fp_promote_t<_Ta>>;
      return emsr::detail::gamma<type>(a);
    }

  // Upper incomplete gamma functions

  /**
   * Return the upper incomplete gamma function @f$ \Gamma(a,x) @f$.
   * The (upper) incomplete gamma function is defined by
   * @f[
   *   \Gamma(a,x) = \int_x^\infty t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename Tp>
    inline emsr::fp_promote_t<_Ta, Tp>
    tgamma(_Ta a, Tp x)
    {
      using type = emsr::fp_promote_t<_Ta, Tp>;
      return emsr::detail::tgamma<type>(a, x);
    }

  // Lower incomplete gamma functions

  /**
   * Return the lower incomplete gamma function @f$ \gamma(a,x) @f$.
   * The lower incomplete gamma function is defined by
   * @f[
   *   \gamma(a,x) = \int_0^x t^{a-1}e^{-t}dt
   * @f]
   */
  template<typename _Ta, typename Tp>
    inline emsr::fp_promote_t<_Ta, Tp>
    tgamma_lower(_Ta a, Tp x)
    {
      using type = emsr::fp_promote_t<_Ta, Tp>;
      return emsr::detail::tgamma_lower<type>(a, x);
    }

  // Digamma or psi functions

  /**
   * Return the digamma or psi function of argument @c x.
   *
   * The the digamma or psi function is defined by
   * @f[
   *    \psi(x) = \frac{d}{dx}log\left(\Gamma(x)\right)
   *            = \frac{\Gamma'(x)}{\Gamma(x)},
   * @f]
   * the logarithmic derivative of the gamma function.
   *
   * @param x The parameter
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    digamma(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::digamma<type>(x);
    }

  // Polygamma or psi functions

  /**
   * Return the polygamma function of argument @c x.
   *
   * The the polygamma or digamma function is defined by
   * @f[
   *    \psi(x) = \frac{d}{dx}log\left(\Gamma(x)\right)
   *            = \frac{\Gamma'(x)}{\Gamma(x)}
   * @f]
   *
   * @param m The order
   * @param x The parameter
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    polygamma(unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::polygamma<type>(m, x);
    }

  // Reciprocal Gamma function.

  /**
   * Return the reciprocal gamma function for real argument.
   *
   * The reciprocal of the Gamma function is what you'd expect:
   * @f[
   *     \Gamma_r(a) = \frac{1}{\Gamma(a)}
   * @f]
   * But unlike the Gamma function this function has no singularities
   * and is exponentially decreasing for increasing argument.
   */
  template<typename _Ta>
    inline emsr::fp_promote_t<_Ta>
    gamma_reciprocal(_Ta a)
    {
      using type = emsr::fp_promote_t<_Ta>;
      return emsr::detail::gamma_reciprocal<type>(a);
    }

  // Scaled lower incomplete gamma

  /**
   * @brief Return the gamma cumulative propability distribution function
   * or the regularized lower incomplete gamma function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename _Ta, typename Tp>
    inline emsr::fp_promote_t<_Ta, Tp>
    gamma_p(_Ta a, Tp x)
    {
      using type = emsr::fp_promote_t<_Ta, Tp>;
      return emsr::detail::gamma_p<type>(a, x);
    }

  // Scaled upper incomplete gamma

  /**
   * @brief Return the gamma complementary cumulative propability distribution
   * (or survival) function or the regularized upper incomplete gamma function.
   *
   * The formula for the gamma probability density function is:
   * @f[
   *    \Gamma(x|\alpha,\beta) = \frac{1}{\beta\Gamma(\alpha)}
   *                             (x/\beta)^{\alpha - 1} e^{-x/\beta}
   * @f]
   */
  template<typename _Ta, typename Tp>
    inline emsr::fp_promote_t<_Ta, Tp>
    gamma_q(_Ta a,Tp  x)
    {
      using type = emsr::fp_promote_t<_Ta, Tp>;
      return emsr::detail::gamma_q<type>(a, x);
    }

  /**
   * Return the harmonic number @f$ H_n @f$.
   *
   * The the harmonic number is defined by
   * @f[
   *    H_n = \sum_{k=1}^{n}\frac{1}{k}
   * @f]
   *
   * @param n The parameter
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    harmonic(unsigned int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::harmonic_number<type>(n);
    }

} // namespace emsr

#endif // SF_GAMMA_H
