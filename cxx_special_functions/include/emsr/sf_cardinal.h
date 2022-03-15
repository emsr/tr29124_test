#ifndef SF_CARDINAL_H
#define SF_CARDINAL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_cardinal.tcc>

namespace emsr
{

  // Sinus cardinal functions

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
