#ifndef SF_AIRY_H
#define SF_AIRY_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_airy.tcc>

namespace emsr
{

  /**
   * Return the Airy function @f$ Ai(x) @f$ of complex argument @c x.
   *
   * The Airy function is defined by:
   * @f[
   *    Ai(x) = \frac{1}{\pi}\int_0^\infty
   *      \cos \left(\frac{t^3}{3} + xt \right)dt
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The complex argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    airy_ai(const std::complex<Tp>& z)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::airy_ai<type>(z);
    }

  /**
   * Return the Airy function @f$ Bi(x) @f$ of complex argument @c x.
   *
   * The Airy function is defined by:
   * @f[
   *    Bi(x) = \frac{1}{\pi}\int_0^\infty \left[ 
   *           \exp \left(-\frac{t^3}{3} + xt \right)
   *          + \sin \left(\frac{t^3}{3} + xt \right) \right] dt
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The complex argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    airy_bi(const std::complex<Tp>& z)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::airy_bi<type>(z);
    }

} // namespace emsr

#endif // SF_AIRY_H
