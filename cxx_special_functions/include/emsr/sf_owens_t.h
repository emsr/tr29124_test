#ifndef SF_OWENS_T_H
#define SF_OWENS_T_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_owens_t.tcc>

namespace emsr
{

  // Owens T functions.

  /**
   * Return the Owens T function @f$ T(h,a) @f$ of shape factor @c h
   * and integration limit @c a.
   *
   * The Owens T function is defined by
   * @f[
   *    T(h,a) = \frac{1}{2\pi}\int_0^a
   *           \frac{\exp\left[-\frac{1}{2}h^2(1+x^2)\right]}{1+x^2} dx
   * @f]
   *
   * @param h The shape factor
   * @param a The integration limit
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    owens_t(Tp h, Tp a)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::owens_t<type>(h, a);
    }

} // namespace emsr

#endif // SF_OWENS_T_H
