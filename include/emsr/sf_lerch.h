#ifndef SF_LERCH_H
#define SF_LERCH_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_lerch.tcc>

namespace emsr
{

  /**
   * Return the Lerch transcendent @f$ \Phi(z,s,a) @f$ for @c float arguments.
   *
   * @see lerch_phi for details.
   */
  inline float
  lerch_phif(float z, float s, float a)
  { return emsr::detail::lerch_phi<float>(z, s, a); }

  /**
   * Return the Lerch transcendent @f$ \Phi(z,s,a) @f$
   * for <tt> long double </tt> arguments.
   *
   * @see lerch_phi for details.
   */
  inline long double
  lerch_phil(long double z, long double s, long double a)
  { return emsr::detail::lerch_phi<long double>(z, s, a); }

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
  template<typename _Tp, typename _Ts, typename _Ta>
    inline emsr::fp_promote_t<_Tp, _Ts, _Ta>
    lerch_phi(_Tp z, _Ts s, _Ta a)
    {
      using type = emsr::fp_promote_t<_Tp, _Ts, _Ta>;
      return emsr::detail::lerch_phi<type>(z, s, a);
    }

} // namespace emsr

#endif // SF_LERCH_H
