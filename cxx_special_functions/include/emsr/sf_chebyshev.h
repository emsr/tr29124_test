#ifndef SF_CHEBYSHEV_H
#define SF_CHEBYSHEV_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_chebyshev.tcc>

namespace emsr
{

  // Chebyshev polynomials of the first kind

  /**
   * Return the Chebyshev polynomial of the first kind @f$ T_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * The Chebyshev polynomial of the first kind is defined by:
   * @f[
   *    T_n(x) = \cos(n \theta)
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    chebyshev_t(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::chebyshev_t<type>(n, x).T_n;
    }

  // Chebyshev polynomials of the second kind

  /**
   * Return the Chebyshev polynomial of the second kind @f$ U_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * The Chebyshev polynomial of the second kind is defined by:
   * @f[
   *    U_n(x) = \frac{\sin \left[(n+1)\theta \right]}{\sin(\theta)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    chebyshev_u(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::chebyshev_u<type>(n, x).U_n;
    }

  // Chebyshev polynomials of the third kind

  /**
   * Return the Chebyshev polynomial of the third kind @f$ V_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * The Chebyshev polynomial of the third kind is defined by:
   * @f[
   *    V_n(x) = \frac{\cos \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\cos \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    chebyshev_v(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::chebyshev_v<type>(n, x).V_n;
    }

  // Chebyshev polynomials of the fourth kind

  /**
   * Return the Chebyshev polynomial of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * The Chebyshev polynomial of the fourth kind is defined by:
   * @f[
   *    W_n(x) = \frac{\sin \left[ \left(n+\frac{1}{2}\right)\theta \right]}
   *                  {\sin \left(\frac{\theta}{2}\right)}
   * @f]
   * where @f$ \theta = \arccos(x) @f$, @f$ -1 <= x <= +1 @f$.
   *
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    chebyshev_w(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::chebyshev_w<type>(n, x).W_n;
    }

} // namespace emsr

#endif // SF_CHEBYSHEV_H
