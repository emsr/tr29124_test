#ifndef SF_FRESNEL_H
#define SF_FRESNEL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_fresnel.tcc>

namespace emsr
{

  // Fresnel sine integral

  inline float
  fresnel_sf(float x)
  { return std::imag(emsr::detail::fresnel<float>(x)); }

  inline long double
  fresnel_sl(long double x)
  { return std::imag(emsr::detail::fresnel<long double>(x)); }

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
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    fresnel_s(_Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return std::imag(emsr::detail::fresnel<type>(x));
    }

  // Fresnel cosine integral

  inline float
  fresnel_cf(float x)
  { return std::real(emsr::detail::fresnel<float>(x)); }

  inline long double
  fresnel_cl(long double x)
  { return std::real(emsr::detail::fresnel<long double>(x)); }

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
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    fresnel_c(_Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return std::real(emsr::detail::fresnel<type>(x));
    }

} // namespace emsr

#endif // SF_FRESNEL_H
