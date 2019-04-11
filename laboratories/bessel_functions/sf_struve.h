
#include "sf_struve.tcc"

namespace __gnu_cxx
{

  // Struve functions (of the first kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline float
  struve_hf(float __nu, float __x)
  { return std::__detail::__struve_h<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{H}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_h for setails.
   */
  inline long double
  struve_hl(long double __nu, long double __x)
  { return std::__detail::__struve_h<long double>(__nu, __x); }

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
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline fp_promote_t<_Tpnu, _Tp>
    struve_h(_Tpnu __nu, _Tp __x)
    {
      using __type = fp_promote_t<_Tpnu, _Tp>;
      return std::__detail::__struve_h<__type>(__nu, __x);
    }

  // Struve functions (of the second kind)

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline float
  struve_kf(float __nu, float __x)
  { return std::__detail::__struve_k<float>(__nu, __x); }

  /**
   * Return the Struve function of the first kind @f$ \boldmath{K}_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see struve_k for setails.
   */
  inline long double
  struve_kl(long double __nu, long double __x)
  { return std::__detail::__struve_k<long double>(__nu, __x); }

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
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline fp_promote_t<_Tpnu, _Tp>
    struve_k(_Tpnu __nu, _Tp __x)
    {
      using __type = fp_promote_t<_Tpnu, _Tp>;
      return std::__detail::__struve_k<__type>(__nu, __x);
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
  struve_lf(float __nu, float __x)
  { return std::__detail::__struve_l<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_l for setails.
   */
  inline long double
  struve_ll(long double __nu, long double __x)
  { return std::__detail::__struve_l<long double>(__nu, __x); }

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
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline fp_promote_t<_Tpnu, _Tp>
    struve_l(_Tpnu __nu, _Tp __x)
    {
      using __type = fp_promote_t<_Tpnu, _Tp>;
      return std::__detail::__struve_l<__type>(__nu, __x);
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
  struve_mf(float __nu, float __x)
  { return std::__detail::__struve_m<float>(__nu, __x); }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_{\nu}(x) @f$ for <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see struve_m for setails.
   */
  inline long double
  struve_ml(long double __nu, long double __x)
  { return std::__detail::__struve_m<long double>(__nu, __x); }

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
   * @tparam _Tpnu The floating-point type of the order @c __nu.
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param  __nu  The order
   * @param  __x   The argument, <tt> __x >= 0 </tt>
   * @throw std::domain_error if <tt> __x < 0 </tt>.
   */
  template<typename _Tpnu, typename _Tp>
    inline fp_promote_t<_Tpnu, _Tp>
    struve_m(_Tpnu __nu, _Tp __x)
    {
      using __type = fp_promote_t<_Tpnu, _Tp>;
      return std::__detail::__struve_m<__type>(__nu, __x);
    }

} // namespace __gnu_cxx
