#ifndef SF_TRIGINT_H
#define SF_TRIGINT_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_trigint.tcc>

namespace emsr
{

  // Sine integrals

  /**
   * Return the sine integral @f$ Si(x) @f$ of @c float argument @c x.
   *
   * @see sinint for details.
   */
  inline float
  sinintf(float x)
  { return emsr::detail::sincosint<float>(x).first; }

  /**
   * Return the sine integral @f$ Si(x) @f$ of <tt>long double</tt>
   * argument @c x.
   *
   * @see sinint for details.
   */
  inline long double
  sinintl(long double x)
  { return emsr::detail::sincosint<long double>(x).first; }

  /**
   * Return the sine integral @f$ Si(x) @f$ of real argument @c x.
   *
   * The sine integral is defined by
   * @f[
   *    Si(x) = \int_0^x \frac{sin(t)}{t}dt
   * @f]
   *
   * @param x The real upper integration limit
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sincosint<type>(x).first;
    }

  // Cosine integrals

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of @c float argument @c x.
   *
   * @see cosint for details.
   */
  inline float
  cosintf(float x)
  { return emsr::detail::sincosint<float>(x).second; }

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of <tt>long double</tt>
   * argument @c x.
   *
   * @see cosint for details.
   */
  inline long double
  cosintl(long double x)
  { return emsr::detail::sincosint<long double>(x).second; }

  /**
   * Return the cosine integral @f$ Ci(x) @f$ of real argument @c x.
   *
   * The cosine integral is defined by
   * @f[
   *    Ci(x) = -\int_x^\infty \frac{cos(t)}{t}dt
   *     = \gamma_E + ln(x) + \int_0^x \frac{cos(t)-1}{t}dt
   * @f]
   *
   * @param x The real upper integration limit
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    cosint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sincosint<type>(x).second;
    }

} // namespace emsr

#endif // SF_TRIGINT_H