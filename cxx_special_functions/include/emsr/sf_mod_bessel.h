#ifndef SF_MOD_BESSEL_H
#define SF_MOD_BESSEL_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_mod_bessel.tcc>

namespace emsr
{

  // Regular modified cylindrical Bessel functions

  /**
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline float
  cyl_bessel_if(float nu, float x)
  { return detail::cyl_bessel_i<float>(nu, x); }

  /**
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_i for setails.
   */
  inline long double
  cyl_bessel_il(long double nu, long double x)
  { return detail::cyl_bessel_i<long double>(nu, x); }

  /**
   * Return the regular modified Bessel function @f$ I_{\nu}(x) @f$
   * for real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = i^{-\nu}J_\nu(ix) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_bessel_i(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return detail::cyl_bessel_i<type>(nu, x);
    }

  // Irregular modified cylindrical Bessel functions

  /**
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline float
  cyl_bessel_kf(float nu, float x)
  { return detail::cyl_bessel_k<float>(nu, x); }

  /**
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_k for setails.
   */
  inline long double
  cyl_bessel_kl(long double nu, long double x)
  { return detail::cyl_bessel_k<long double>(nu, x); }

  /**
   * Return the irregular modified Bessel function @f$ K_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @c x.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_bessel_k(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return detail::cyl_bessel_k<type>(nu, x);
    }

  // Scaled modified Bessel functions

  /**
   * Return the scaled regular modified Bessel function
   * @f$ e^{-x} I_{\nu}(x) @f$
   * for real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The regular modified cylindrical Bessel function is:
   * @f[
   *  I_{\nu}(x) = i^{-\nu}J_\nu(ix) = \sum_{k=0}^{\infty}
   * 		\frac{(x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_bessel_i_scaled(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_bessel_i_scaled<type>(nu, x);
    }

  /**
   * Return the scaled irregular modified Bessel function
   * @f$ e^{x} K_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x @f$.
   *
   * The irregular modified Bessel function is defined by:
   * @f[
   * 	K_{\nu}(x) = \frac{\pi}{2}
   * 		     \frac{I_{-\nu}(x) - I_{\nu}(x)}{\sin \nu\pi}
   * @f]
   * where for integral @f$ \nu = n @f$ a limit is taken:
   * @f$ lim_{\nu \to n} @f$.
   * For negative argument we have simply:
   * @f[
   * 	K_{-\nu}(x) = K_{\nu}(x)
   * @f]
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_bessel_k_scaled(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_bessel_k_scaled<type>(nu, x);
    }

  // Modified spherical Bessel functions of the first kind

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_i for details.
   */
  inline float
  sph_bessel_if(unsigned int n, float x)
  { return emsr::detail::sph_bessel_ik<float>(n, x).i_value; }

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_i for details.
   */
  inline long double
  sph_bessel_il(unsigned int n, long double x)
  { return emsr::detail::sph_bessel_ik<long double>(n, x).i_value; }

  /**
   * Return the regular modified spherical Bessel function @f$ i_n(x) @f$
   * of nonnegative order n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  i_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} I_{n+1/2}(x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sph_bessel_i(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_bessel_ik<type>(n, x).i_value;
    }

  // Modified spherical Bessel functions of the second kind

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_k for more details.
   */
  inline float
  sph_bessel_kf(unsigned int n, float x)
  { return emsr::detail::sph_bessel_ik<float>(n, x).k_value; }

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel_k for more details.
   */
  inline long double
  sph_bessel_kl(unsigned int n, long double x)
  { return emsr::detail::sph_bessel_ik<long double>(n, x).k_value; }

  /**
   * Return the irregular modified spherical Bessel function @f$ k_n(x) @f$
   * of nonnegative order n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  k_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} K_{n+1/2}(x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sph_bessel_k(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_bessel_ik<type>(n, x).k_value;
    }

  // Airy functions of the first kind

  /**
   * Return the Airy function @f$ Ai(x) @f$ for @c float argument @c x.
   *
   * @see airy_ai for details.
   */
  inline float
  airy_aif(float x)
  { return emsr::detail::airy<float>(x).Ai_value; }

  /**
   * Return the Airy function @f$ Ai(x) @f$ for <tt>long double</tt>
   * argument @c x.
   *
   * @see airy_ai for details.
   */
  inline long double
  airy_ail(long double x)
  { return emsr::detail::airy<long double>(x).Ai_value; }

  /**
   * Return the Airy function @f$ Ai(x) @f$ of real argument @c x.
   *
   * The Airy function is defined by:
   * @f[
   *    Ai(x) = \frac{1}{\pi}\int_0^\infty
   *      \cos \left(\frac{t^3}{3} + xt \right)dt
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    airy_ai(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::airy<type>(x).Ai_value;
    }

  // Airy functions of the second kind

  /**
   * Return the Airy function @f$ Bi(x) @f$ for @c float argument @c x.
   *
   * @see airy_bi for details.
   */
  inline float
  airy_bif(float x)
  { return emsr::detail::airy<float>(x).Bi_value; }

  /**
   * Return the Airy function @f$ Bi(x) @f$ for <tt>long double</tt>
   * argument @c x.
   *
   * @see airy_bi for details.
   */
  inline long double
  airy_bil(long double x)
  { return emsr::detail::airy<long double>(x).Bi_value; }

  /**
   * Return the Airy function @f$ Bi(x) @f$ of real argument @c x.
   *
   * The Airy function is defined by:
   * @f[
   *    Bi(x) = \frac{1}{\pi}\int_0^\infty \left[
   *           \exp \left(-\frac{t^3}{3} + xt \right)
   *          + \sin \left(\frac{t^3}{3} + xt \right) \right] dt
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    airy_bi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::airy<type>(x).Bi_value;
    }

} // namespace emsr

#endif // SF_MOD_BESSEL_H