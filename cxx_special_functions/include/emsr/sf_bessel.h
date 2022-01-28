#ifndef SF_BESSEL_H
#define SF_BESSEL_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_bessel.tcc>

namespace emsr
{

  // Cylindrical Bessel functions (of the first kind)

  /**
   * Return the Bessel function of the first kind @f$ J_{\nu}(x) @f$
   * for @c float order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline float
  cyl_bessel_jf(float nu, float x)
  { return detail::cyl_bessel_j<float>(nu, x); }

  /**
   * Return the Bessel function of the first kind @f$ J_{\nu}(x) @f$
   * for <tt>long double</tt> order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * @see cyl_bessel_j for setails.
   */
  inline long double
  cyl_bessel_jl(long double nu, long double x)
  { return detail::cyl_bessel_j<long double>(nu, x); }

  /**
   * Return the Bessel function @f$ J_{\nu}(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The cylindrical Bessel function is:
   * @f[
   *    J_{\nu}(x) = \sum_{k=0}^{\infty}
   *              \frac{(-1)^k (x/2)^{\nu + 2k}}{k!\Gamma(\nu+k+1)}
   * @f]
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename _Tp>
    inline emsr::fp_promote_t<Tpnu, _Tp>
    cyl_bessel_j(Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, _Tp>;
      return detail::cyl_bessel_j<type>(nu, x);
    }

  // Cylindrical Neumann functions

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of @c float order @f$ \nu @f$ and argument @c x.
   *
   * @see cyl_neumann for setails.
   */
  inline float
  cyl_neumannf(float nu, float x)
  { return detail::cyl_neumann_n<float>(nu, x); }

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of <tt>long double</tt> order @f$ \nu @f$ and argument @c x.
   *
   * @see cyl_neumann for setails.
   */
  inline long double
  cyl_neumannl(long double nu, long double x)
  { return detail::cyl_neumann_n<long double>(nu, x); }

  /**
   * Return the Neumann function @f$ N_{\nu}(x) @f$
   * of real order @f$ \nu @f$ and argument @f$ x >= 0 @f$.
   *
   * The Neumann function is defined by:
   * @f[
   *    N_{\nu}(x) = \frac{J_{\nu}(x) \cos \nu\pi - J_{-\nu}(x)}
   *                      {\sin \nu\pi}
   * @f]
   * where @f$ x >= 0 @f$ and for integral order @f$ \nu = n @f$
   * a limit is taken: @f$ lim_{\nu \to n} @f$.
   *
   * @tparam Tpnu The floating-point type of the order @c nu.
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  nu  The order.
   * @param  x   The argument, <tt> x >= 0 </tt>.
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename _Tp>
    inline emsr::fp_promote_t<Tpnu, _Tp>
    cyl_neumann(Tpnu nu, _Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, _Tp>;
      return detail::cyl_neumann_n<type>(nu, x);
    }

  // Spherical Bessel functions

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel for more details.
   */
  inline float
  sph_besself(unsigned int n, float x)
  { return detail::sph_bessel<float>(n, x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_bessel for more details.
   */
  inline long double
  sph_bessell(unsigned int n, long double x)
  { return detail::sph_bessel<long double>(n, x); }

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    sph_bessel(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return detail::sph_bessel<type>(n, x);
    }

  // Spherical Neumann functions

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_neumann for details.
   */
  inline float
  sph_neumannf(unsigned int n, float x)
  { return detail::sph_neumann<float>(n, x); }

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and <tt>long double</tt> @f$ x >= 0 @f$.
   *
   * @see sph_neumann for details.
   */
  inline long double
  sph_neumannl(unsigned int n, long double x)
  { return detail::sph_neumann<long double>(n, x); }

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Neumann function is defined by
   * @f[
   *    n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    sph_neumann(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return detail::sph_neumann<type>(n, x);
    }

  // Cylindrical Hankel functions of the first kind

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_1 for details.
   */
  inline std::complex<float>
  cyl_hankel_1f(float nu, float z)
  { return emsr::detail::cyl_hankel_1<float>(nu, z); }

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_1 for details.
   */
  inline std::complex<long double>
  cyl_hankel_1l(long double nu, long double z)
  { return emsr::detail::cyl_hankel_1<long double>(nu, z); }

  /**
   * Return the cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_n(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The spherical Hankel function of the first kind is defined by:
   * @f[
   *    H^{(1)}_\nu(x) = J_\nu(x) + iN_\nu(x)
   * @f]
   * where @f$ J_\nu(x) @c and N_\nu(x) @f$ are the cylindrical Bessel
   * and Neumann functions respectively (see cyl_bessel and cyl_neumann).
   *
   * @tparam _Tp The real type of the argument
   * @param nu The real order
   * @param z The real argument
   */
  template<typename Tpnu, typename _Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, _Tp>>
    cyl_hankel_1(Tpnu nu, _Tp z)
    {
      using type = emsr::fp_promote_t<Tpnu, _Tp>;
      return emsr::detail::cyl_hankel_1<type>(nu, z);
    }

  // Cylindrical Hankel functions of the second kind

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of @c float order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_2 for details.
   */
  inline std::complex<float>
  cyl_hankel_2f(float nu, float z)
  { return emsr::detail::cyl_hankel_2<float>(nu, z); }

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>long double</tt> order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * @see cyl_hankel_2 for details.
   */
  inline std::complex<long double>
  cyl_hankel_2l(long double nu, long double z)
  { return emsr::detail::cyl_hankel_2<long double>(nu, z); }

  /**
   * Return the cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_n(x) @f$ of real order @f$ \nu @f$
   * and argument @f$ x >= 0 @f$.
   *
   * The cylindrical Hankel function of the second kind is defined by:
   * @f[
   *    H^{(2)}_\nu(x) = J_\nu(x) - iN_\nu(x)
   * @f]
   * where @f$ J_\nu(x) @c and N_\nu(x) @f$ are the cylindrical Bessel
   * and Neumann functions respectively (see cyl_bessel and cyl_neumann).
   *
   * @tparam _Tp The real type of the argument
   * @param nu The real order
   * @param z The real argument
   */
  template<typename Tpnu, typename _Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, _Tp>>
    cyl_hankel_2(Tpnu nu, _Tp z)
    {
      using type = emsr::fp_promote_t<Tpnu, _Tp>;
      return emsr::detail::cyl_hankel_2<type>(nu, z);
    }

  // Spherical Hankel functions of the first kind

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_1 for details.
   */
  inline std::complex<float>
  sph_hankel_1f(unsigned int n, float z)
  { return emsr::detail::sph_hankel_1<float>(n, z); }

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order n and @c <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_1 for details.
   */
  inline std::complex<long double>
  sph_hankel_1l(unsigned int n, long double z)
  { return emsr::detail::sph_hankel_1<long double>(n, z); }

  /**
   * Return the spherical Hankel function of the first kind @f$ h^{(1)}_n(x) @f$
   * of nonnegative order @c n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Hankel function of the first kind is defined by:
   * @f[
   *    h^{(1)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(1)}_{n+1/2}(x)
   * @f]
   * or in terms of the cylindrical Bessel and Neumann functions by:
   * @f[
   *    h^{(1)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2}
   *       \left[ J_{n+1/2}(x) + iN_{n+1/2}(x) \right]
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param n The non-negative order
   * @param z The real argument
   */
  template<typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tp>>
    sph_hankel_1(unsigned int n, _Tp z)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::sph_hankel_1<type>(n, z);
    }

  // Spherical Hankel functions of the second kind

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order n and @c float argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_2 for details.
   */
  inline std::complex<float>
  sph_hankel_2f(unsigned int n, float z)
  { return emsr::detail::sph_hankel_2<float>(n, z); }

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order n and @c <tt>long double</tt> argument @f$ x >= 0 @f$.
   *
   * @see sph_hankel_2 for details.
   */
  inline std::complex<long double>
  sph_hankel_2l(unsigned int n, long double z)
  { return emsr::detail::sph_hankel_2<long double>(n, z); }

  /**
   * Return the spherical Hankel function of the second kind @f$ h^{(2)}_n(x)@f$
   * of nonnegative order @c n and real argument @f$ x >= 0 @f$.
   *
   * The spherical Hankel function of the second kind is defined by:
   * @f[
   *    h^{(2)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(2)}_{n+1/2}(x)
   * @f]
   * or in terms of the cylindrical Bessel and Neumann functions by:
   * @f[
   *    h^{(2)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2}
   *       \left[ J_{n+1/2}(x) - iN_{n+1/2}(x) \right]
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param n The non-negative order
   * @param z The real argument
   */
  template<typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tp>>
    sph_hankel_2(unsigned int n, _Tp z)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::sph_hankel_2<type>(n, z);
    }

} // namespace emsr

#endif // SF_BESSEL_H
