#ifndef SF_FRESNEL_H
#define SF_FRESNEL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_fresnel.tcc>

namespace emsr
{

  // Fresnel sine integral

  /**
   * Return the Fresnel sine integral of argument @c x.
   *
   * The Fresnel sine integral is defined by
   * @f[
   *    S(x) = \int_0^x \sin(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    fresnel_s(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return std::imag(emsr::detail::fresnel<type>(x));
    }

  // Fresnel cosine integral

  /**
   * Return the Fresnel cosine integral of argument @c x.
   *
   * The Fresnel cosine integral is defined by
   * @f[
   *    C(x) = \int_0^x \cos(\frac{\pi}{2}t^2) dt
   * @f]
   *
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    fresnel_c(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return std::real(emsr::detail::fresnel<type>(x));
    }

} // namespace emsr

#endif // SF_FRESNEL_H
