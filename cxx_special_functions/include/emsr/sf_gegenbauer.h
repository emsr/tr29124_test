#ifndef SF_GEGENBAUER_H
#define SF_GEGENBAUER_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_gegenbauer.tcc>

namespace emsr
{

  // Gegenbauer polynomials

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{(\lambda)}(x) @f$ of degree @c n
   * and @c float order @f$ \lambda > -1/2, \lambda \neq 0 @f$
   * and argument @c x.
   *
   * @see gegenbauer for details.
   */
  inline float
  gegenbauerf(unsigned int n, float lambda, float x)
  { return emsr::detail::gegenbauer_recur<float>(n, lambda, x).C_n; }

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{\lambda}(x) @f$ of degree @c n
   * and <tt>long double</tt> order @f$ \lambda > -1/2, \lambda \neq 0 @f$
   * and argument @c x.
   *
   * @see gegenbauer for details.
   */
  inline long double
  gegenbauerl(unsigned int n, long double lambda, long double x)
  {
    return emsr::detail::gegenbauer_recur<long double>(n, lambda, x).C_n;
  }

  /**
   * Return the Gegenbauer polynomial @f$ C_n^{\lambda}(x) @f$ of degree @c n
   * and real order @f$ \lambda > -1/2, \lambda \neq 0 @f$ and argument
   * @c x.
   *
   * The Gegenbauer polynomial is generated by a three-term recursion relation:
   * @f[
   *    C_n^{\lambda}(x) = \frac{1}{n}\left[ 2x(n+\lambda-1)C_{n-1}^{\lambda}(x)
   *                   - (n+2\lambda-2)C_{n-2}^{\lambda}(x) \right]
   * @f]
   * and @f$ C_0^{\lambda}(x) = 1 @f$, @f$ C_1^{\lambda}(x) = 2\lambda x @f$.
   *
   * @tparam Tlam The real type of the order
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral degree
   * @param lambda The real order
   * @param x The real argument
   */
  template<typename Tlam, typename Tp>
    inline typename emsr::fp_promote_t<Tlam, Tp>
    gegenbauer(unsigned int n, Tlam lambda, Tp x)
    {
      using type = emsr::fp_promote_t<Tlam, Tp>;
      return emsr::detail::gegenbauer_recur<type>(n, lambda, x).C_n;
    }

} // namespace emsr

#endif // SF_GEGENBAUER_H