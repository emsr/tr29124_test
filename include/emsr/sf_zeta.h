#ifndef SF_ZETA_H
#define SF_ZETA_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_zeta.tcc>

namespace emsr
{

  // Riemann zeta functions

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for @c float argument @c s.
   *
   * @see riemann_zeta for more details.
   */
  inline float
  riemann_zetaf(float s)
  { return emsr::detail::riemann_zeta<float>(s); }

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for <tt>long double</tt> argument @c s.
   *
   * @see riemann_zeta for more details.
   */
  inline long double
  riemann_zetal(long double s)
  { return emsr::detail::riemann_zeta<long double>(s); }

  /**
   * Return the Riemann zeta function @f$ \zeta(s) @f$
   * for real argument @c s.
   *
   * The Riemann zeta function is defined by:
   * @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} k^{-s} \mbox{ for } s > 1
   * @f]
   * and
   * @f[
   * 	\zeta(s) = \frac{1}{1-2^{1-s}}\sum_{k=1}^{\infty}(-1)^{k-1}k^{-s}
   *              \mbox{ for } 0 <= s < 1
   * @f]
   * For s < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = 2^s \pi^{s-1} \sin(\frac{\pi s}{2}) \Gamma(1-s) \zeta(1-s)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c s.
   * @param s The order <tt> s != 1 </tt>
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    riemann_zeta(Tp s)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::riemann_zeta<type>(s);
    }

  // Hurwitz zeta functions

  /**
   * Return the Hurwitz zeta function of @c float argument @c s,
   * and parameter @c a.
   *
   * @see hurwitz_zeta for details.
   */
  inline float
  hurwitz_zetaf(float s, float a)
  { return emsr::detail::hurwitz_zeta<float>(s, a); }

  /**
   * Return the Hurwitz zeta function of <tt>long double</tt>
   * argument @c s, and parameter @c a.
   *
   * @see hurwitz_zeta for details.
   */
  inline long double
  hurwitz_zetal(long double s, long double a)
  { return emsr::detail::hurwitz_zeta<long double>(s, a); }

  /**
   * Return the Hurwitz zeta function of real argument @c s,
   * and parameter @c a.
   *
   * The the Hurwitz zeta function is defined by
   * @f[
   *    \zeta(s, a) = \sum_{n=0}^{\infty}\frac{1}{(a + n)^s}
   * @f]
   *
   * @param s The order.
   * @param a The parameter.
   */
  template<typename Tp, typename _Up>
    inline emsr::fp_promote_t<Tp, _Up>
    hurwitz_zeta(Tp s, _Up a)
    {
      using type = emsr::fp_promote_t<Tp, _Up>;
      return emsr::detail::hurwitz_zeta<type>(s, a);
    }

  // Dilogarithm functions

  /**
   * Return the dilogarithm function @f$ Li_2(z) @f$
   * for @c float argument.
   *
   * @see dilog for details.
   */
  inline float
  dilogf(float x)
  { return emsr::detail::dilog<float>(x); }

  /**
   * Return the dilogarithm function @f$ Li_2(z) @f$
   * for <tt>long double</tt> argument.
   *
   * @see dilog for details.
   */
  inline long double
  dilogl(long double x)
  { return emsr::detail::dilog<long double>(x); }

  /**
   * Return the dilogarithm function @f$ Li_2(z) @f$
   * for real argument.
   *
   * The dilogarithm is defined by:
   * @f[
   *    Li_2(x) = \sum_{k=1}^{\infty}\frac{x^k}{k^2}
   * @f]
   *
   * @param x The argument.
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    dilog(_Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::dilog<type>(x);
    }

  // Debye functions.

  /**
   * Return the Debye function @f$ D_n(x) @f$
   * of positive order @c n and @c float argument @c x.
   *
   * @see debye for details.
   */
  inline float
  debyef(unsigned int n, float x)
  { return emsr::detail::debye<float>(n, x); }

  /**
   * Return the Debye function @f$ D_n(x) @f$
   * of positive order @c n and real argument @c x.
   *
   * @see debye for details.
   */
  inline long double
  debyel(unsigned int n, long double x)
  { return emsr::detail::debye<long double>(n, x); }

  /**
   * Return the Debye function @f$ D_n(x) @f$
   * of positive order @c n and real argument @c x.
   *
   * The Debye function is defined by:
   * @f[
   *    D_n(x) = \frac{n}{x^n}\int_{0}^{x}\frac{t^n}{e^t-1}dt
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param n The positive integral order
   * @param x The real argument @f$ x >= 0 @f$
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    debye(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::debye<type>(n, x);
    }

} // namespace emsr

#endif // SF_ZETA_H
