#ifndef SF_LAGUERRE_H
#define SF_LAGUERRE_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_laguerre.tcc>

namespace emsr
{

  // Associated Laguerre polynomials

  /**
   * Return the associated Laguerre polynomial @f$ L_n^{(m)}(x) @f$
   * of order @c n, degree @c m, and @c float argument @c x.
   *
   * @see assoc_laguerre for more details.
   */
  inline float
  assoc_laguerref(unsigned int n, unsigned int m, float x)
  { return emsr::detail::assoc_laguerre<float>(n, m, x); }

  /**
   * Return the associated Laguerre polynomial @f$ L_n^{(m)}(x) @f$
   * of order @c n, degree @c m and <tt>long double</tt>
   * argument @c x.
   *
   * @see assoc_laguerre for more details.
   */
  inline long double
  assoc_laguerrel(unsigned int n, unsigned int m, long double x)
  { return emsr::detail::assoc_laguerre<long double>(n, m, x); }

  /**
   * Return the associated Laguerre polynomial @f$ L_n^{(m)}(x) @f$
   * of nonnegative degree @c n, nonnegative order @c m
   * and real argument @c x.
   *
   * The associated Laguerre function of real order @f$ \alpha @f$,
   * @f$ L_n^{(\alpha)}(x) @f$, is defined by
   * @f[
   * 	 L_n^{(\alpha)}(x) = \frac{(\alpha + 1)_n}{n!}
   * 			 {}_1F_1(-n; \alpha + 1; x)
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
   * and @f$ x >= 0 @f$.
   * @see laguerre for details of the Laguerre function of degree @c n
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param n The degree of the Laguerre function, <tt>n >= 0</tt>.
   * @param m The order of the Laguerre function, <tt>m >= 0</tt>.
   * @param x The argument of the Laguerre function, <tt>x >= 0</tt>.
   * @throw std::domain_error if <tt>x < 0</tt>.
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    assoc_laguerre(unsigned int n, unsigned int m, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::assoc_laguerre<type>(n, m, x);
    }

  /**
   * Return the associated Laguerre polynomial @f$ L_n^{(\alpha)}(x) @f$
   * of nonnegative degree @c n, order @f$ \alpha @f$
   * and real argument @c x.
   *
   * @tparam _Talpha The (signed integer or floating-point) type
   *                 of the degree @c alpha1.
   * @tparam _Tp The floating-point type of the argument @c x.
   */
  template<typename _Talpha, typename _Tp>
    inline emsr::fp_promote_t<_Talpha, _Tp>
    assoc_laguerre(unsigned int n, _Talpha alpha1, _Tp x)
    {
      using type = emsr::fp_promote_t<_Talpha, _Tp>;
      return emsr::detail::assoc_laguerre<type>(n, alpha1, x);
    }

  // Laguerre polynomials

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$ of nonnegative degree @c n
   * and @c float argument  @f$ x >= 0 @f$.
   *
   * @see laguerre for more details.
   */
  inline float
  laguerref(unsigned int n, float x)
  { return emsr::detail::laguerre<float>(n, x); }

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$ of nonnegative degree @c n
   * and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see laguerre for more details.
   */
  inline long double
  laguerrel(unsigned int n, long double x)
  { return emsr::detail::laguerre<long double>(n, x); }

  /**
   * Returns the Laguerre polynomial @f$ L_n(x) @f$
   * of nonnegative degree @c n and real argument @f$ x >= 0 @f$.
   *
   * The Laguerre polynomial is defined by:
   * @f[
   * 	 L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} (x^ne^{-x})
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param n The nonnegative order
   * @param x The argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    laguerre(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::laguerre<type>(n, x);
    }

} // namespace emsr

#endif // SF_LAGUERRE_H
