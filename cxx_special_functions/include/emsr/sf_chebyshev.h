#ifndef SF_CHEBYSHEV_H
#define SF_CHEBYSHEV_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_chebyshev.tcc>

namespace emsr
{

  // Chebyshev polynomials of the first kind

  /**
   * Return the Chebyshev polynomials of the first kind @f$ T_n(x) @f$
   * of non-negative order @c n and @c float argument @c x.
   *
   * @see chebyshev_t for details.
   */
  inline float
  chebyshev_tf(unsigned int n, float x)
  { return emsr::detail::chebyshev_t<float>(n, x).T_n; }

  /**
   * Return the Chebyshev polynomials of the first kind @f$ T_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * @see chebyshev_t for details.
   */
  inline long double
  chebyshev_tl(unsigned int n, long double x)
  { return emsr::detail::chebyshev_t<long double>(n, x).T_n; }

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
   * @tparam _Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    chebyshev_t(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::chebyshev_t<type>(n, x).T_n;
    }

  // Chebyshev polynomials of the second kind

  /**
   * Return the Chebyshev polynomials of the second kind @f$ U_n(x) @f$
   * of non-negative order @c n and @c float argument @c x.
   *
   * @see chebyshev_u for details.
   */
  inline float
  chebyshev_uf(unsigned int n, float x)
  { return emsr::detail::chebyshev_u<float>(n, x).U_n; }

  /**
   * Return the Chebyshev polynomials of the second kind  @f$ U_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * @see chebyshev_u for details.
   */
  inline long double
  chebyshev_ul(unsigned int n, long double x)
  { return emsr::detail::chebyshev_u<long double>(n, x).U_n; }

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
   * @tparam _Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    chebyshev_u(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::chebyshev_u<type>(n, x).U_n;
    }

  // Chebyshev polynomials of the third kind

  /**
   * Return the Chebyshev polynomials of the third kind @f$ V_n(x) @f$
   * of non-negative order @c n and @c float argument @c x.
   *
   * @see chebyshev_v for details.
   */
  inline float
  chebyshev_vf(unsigned int n, float x)
  { return emsr::detail::chebyshev_v<float>(n, x).V_n; }

  /**
   * Return the Chebyshev polynomials of the third kind @f$ V_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * @see chebyshev_v for details.
   */
  inline long double
  chebyshev_vl(unsigned int n, long double x)
  { return emsr::detail::chebyshev_v<long double>(n, x).V_n; }

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
   * @tparam _Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    chebyshev_v(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::chebyshev_v<type>(n, x).V_n;
    }

  // Chebyshev polynomials of the fourth kind

  /**
   * Return the Chebyshev polynomials of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @c n and @c float argument @c x.
   *
   * @see chebyshev_w for details.
   */
  inline float
  chebyshev_wf(unsigned int n, float x)
  { return emsr::detail::chebyshev_w<float>(n, x).W_n; }

  /**
   * Return the Chebyshev polynomials of the fourth kind @f$ W_n(x) @f$
   * of non-negative order @c n and real argument @c x.
   *
   * @see chebyshev_w for details.
   */
  inline long double
  chebyshev_wl(unsigned int n, long double x)
  { return emsr::detail::chebyshev_w<long double>(n, x).W_n; }

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
   * @tparam _Tp The real type of the argument
   * @param n The non-negative integral order
   * @param x The real argument @f$ -1 <= x <= +1 @f$
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    chebyshev_w(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::chebyshev_w<type>(n, x).W_n;
    }

} // namespace emsr

#endif // SF_CHEBYSHEV_H
