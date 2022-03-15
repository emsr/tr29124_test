
#include <sf_struve.tcc>

namespace emsr
{

  // Struve functions (of the first kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline float
  struve_hf(float nu, float x)
  { return emsr::detail::struve_h<float>(nu, x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline long double
  struve_hl(long double nu, long double x)
  { return emsr::detail::struve_h<long double>(nu, x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Struve function is:
   * @f[
   *    \boldmath{H}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline emsr::fp_promote_t<_Tpnu, _Tp>
    struve_h(_Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::struve_h<type>(nu, x);
    }

  // Struve functions (of the second kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline float
  struve_kf(float nu, float x)
  { return emsr::detail::struve_k<float>(nu, x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline long double
  struve_kl(long double nu, long double x)
  { return emsr::detail::struve_k<long double>(nu, x); }

  /**
   * Return the Struve function of the second kind @f$ \boldmath{K}_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Struve function of the second kind is:
   * @f[
   *    \boldmath{K}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline emsr::fp_promote_t<_Tpnu, _Tp>
    struve_k(_Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::struve_k<type>(nu, x);
    }

  // Modified Struve functions (of the first kind)

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ for @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline float
  struve_lf(float nu, float x)
  { return emsr::detail::struve_l<float>(nu, x); }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline long double
  struve_ll(long double nu, long double x)
  { return emsr::detail::struve_l<long double>(nu, x); }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    \boldmath{L}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline emsr::fp_promote_t<_Tpnu, _Tp>
    struve_l(_Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::struve_l<type>(nu, x);
    }

  // Modified Struve functions of the second kind

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ for @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline float
  struve_mf(float nu, float x)
  { return emsr::detail::struve_m<float>(nu, x); }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline long double
  struve_ml(long double nu, long double x)
  { return emsr::detail::struve_m<long double>(nu, x); }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The modified Struve function is:
   * @f[
   *    \boldmath{M}_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam _Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline emsr::fp_promote_t<_Tpnu, _Tp>
    struve_m(_Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::struve_m<type>(nu, x);
    }

} // namespace emsr
