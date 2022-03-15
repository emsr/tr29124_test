#ifndef SF_BESSEL_H
#define SF_BESSEL_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_bessel.tcc>

namespace emsr
{

  // Cylindrical Bessel functions (of the first kind)

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
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order
   * @param  x   The argument, <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_bessel_j(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return detail::cyl_bessel_j<type>(nu, x);
    }

  // Cylindrical Neumann functions

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
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  nu  The order.
   * @param  x   The argument, <tt> x >= 0 </tt>.
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    cyl_neumann(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return detail::cyl_neumann_n<type>(nu, x);
    }

  // Spherical Bessel functions

  /**
   * Return the spherical Bessel function @f$ j_n(x) @f$ of nonnegative order n
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Bessel function is defined by:
   * @f[
   *  j_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} J_{n+1/2}(x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sph_bessel(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return detail::sph_bessel<type>(n, x);
    }

  // Spherical Neumann functions

  /**
   * Return the spherical Neumann function of integral order @f$ n >= 0 @f$
   * and real argument @f$ x >= 0 @f$.
   *
   * The spherical Neumann function is defined by
   * @f[
   *    n_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} N_{n+1/2}(x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  n  The integral order <tt> n >= 0 </tt>
   * @param  x  The real argument <tt> x >= 0 </tt>
   * @throw std::domain_error if <tt> x < 0 </tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sph_neumann(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return detail::sph_neumann<type>(n, x);
    }

  // Cylindrical Hankel functions of the first kind

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
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param z The real argument
   */
  template<typename Tpnu, typename Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, Tp>>
    cyl_hankel_1(Tpnu nu, Tp z)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_hankel_1<type>(nu, z);
    }

  // Cylindrical Hankel functions of the second kind

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
   * @tparam Tp The real type of the argument
   * @param nu The real order
   * @param z The real argument
   */
  template<typename Tpnu, typename Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, Tp>>
    cyl_hankel_2(Tpnu nu, Tp z)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_hankel_2<type>(nu, z);
    }

  // Spherical Hankel functions of the first kind

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
   * @tparam Tp The real type of the argument
   * @param n The non-negative order
   * @param z The real argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    sph_hankel_1(unsigned int n, Tp z)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_hankel_1<type>(n, z);
    }

  // Spherical Hankel functions of the second kind

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
   * @tparam Tp The real type of the argument
   * @param n The non-negative order
   * @param z The real argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    sph_hankel_2(unsigned int n, Tp z)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_hankel_2<type>(n, z);
    }

} // namespace emsr

#endif // SF_BESSEL_H
