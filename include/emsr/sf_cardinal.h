#ifndef SF_CARDINAL_H
#define SF_CARDINAL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_cardinal.tcc>

namespace emsr
{

  // Sinus cardinal functions

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for @c float argument @c x.
   *
   * @see sinc_pi for details.
   */
  inline float
  sincf(float x)
  { return emsr::detail::sinc<float>(x); }

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see sinc_pi for details.
   */
  inline long double
  sincl(long double x)
  { return emsr::detail::sinc<long double>(x); }

  /**
   * Return the sinus cardinal function @f$ sinc_\pi(x) @f$
   * for real argument @c x.
   * The sinus cardinal function is defined by:
   * @f[
   *    sinc(x) = \frac{sin(x)}{x}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinc(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinc<type>(x);
    }

  // Normalized sinus cardinal functions

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for @c float argument @c x.
   *
   * @see sinc for details.
   */
  inline float
  sinc_pif(float x)
  { return emsr::detail::sinc_pi<float>(x); }

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see sinc for details.
   */
  inline long double
  sinc_pil(long double x)
  { return emsr::detail::sinc_pi<long double>(x); }

  /**
   * Return the reperiodized sinus cardinal function @f$ sinc(x) @f$
   * for real argument @c x.
   * The normalized sinus cardinal function is defined by:
   * @f[
   *    sinc_\pi(x) = \frac{sin(\pi x)}{\pi x}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinc_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinc_pi<type>(x);
    }

  // Unnormalized hyperbolic sinus cardinal functions

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for @c float argument @c x.
   *
   * @see sinhc_pi for details.
   */
  inline float
  sinhc_pif(float x)
  { return emsr::detail::sinhc_pi<float>(x); }

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see sinhc_pi for details.
   */
  inline long double
  sinhc_pil(long double x)
  { return emsr::detail::sinhc_pi<long double>(x); }

  /**
   * Return the hyperbolic sinus cardinal function @f$ sinhc_\pi(x) @f$
   * for real argument @c x.
   * The sinus cardinal function is defined by:
   * @f[
   *    sinhc_\pi(x) = \frac{\sinh(x)}{x}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinhc_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinhc_pi<type>(x);
    }

  // Normalized hyperbolic sinus cardinal functions

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for @c float argument @c x.
   *
   * @see sinhc for details.
   */
  inline float
  sinhcf(float x)
  { return emsr::detail::sinhc<float>(x); }

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see sinhc for details.
   */
  inline long double
  sinhcl(long double x)
  { return emsr::detail::sinhc<long double>(x); }

  /**
   * Return the normalized hyperbolic sinus cardinal function @f$ sinhc(x) @f$
   * for real argument @c x.
   * The normalized hyperbolic sinus cardinal function is defined by:
   * @f[
   *    sinhc(x) = \frac{\sinh(\pi x)}{\pi x}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinhc(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinhc<type>(x);
    }
} // namespace emsr

#endif // SF_CARDINAL_H
