#ifndef SF_LERCH_H
#define SF_LERCH_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_lerch.tcc>

namespace emsr
{

  /**
   * Return the Lerch transcendent @f$ \Phi(z,s,a) @f$.
   *
   * The series is:
   * @f[
   *   \Phi(z,s,a) = \sum_{k=0}^{\infty}\frac{z^k}{(k+a)^s}
   * @f]
   *
   * @param z The argument.
   * @param s The order @f$ s != 1 @f$.
   * @param a The scale parameter @f$ a > -1 @f$.
   */
  template<typename Tp, typename Ts, typename Ta>
    inline emsr::fp_promote_t<Tp, Ts, Ta>
    lerch_phi(Tp z, Ts s, Ta a)
    {
      using type = emsr::fp_promote_t<Tp, Ts, Ta>;
      return emsr::detail::lerch_phi<type>(z, s, a);
    }

} // namespace emsr

#endif // SF_LERCH_H
