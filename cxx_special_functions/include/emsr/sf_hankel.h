#ifndef SF_HANKEL_H
#define SF_HANKEL_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_hankel.tcc>

namespace emsr
{

  // Cylindrical Hankel functions of the first kind.

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of complex order @f$ \nu @f$
   * and argument @c x.
   *
   * The cylindrical Hankel function of the first kind is defined by
   * @f[
   *    H^{(1)}_\nu(x) = J_\nu(x) + i N_\nu(x)
   * @f]
   *
   * @tparam Tpnu The complex type of the order
   * @tparam Tp The complex type of the argument
   * @param nu The complex order
   * @param x The complex argument
   */
  template<typename Tpnu, typename Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, Tp>>
    cyl_hankel_1(std::complex<Tpnu> nu, std::complex<Tp> x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_hankel_1<type>(nu, x);
    }

  // Cylindrical Hankel functions of the second kind.

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of complex order @f$ \nu @f$
   * and argument @c x.
   *
   * The cylindrical Hankel function of the second kind is defined by
   * @f[
   *    H^{(2)}_\nu(x) = J_\nu(x) - i N_\nu(x)
   * @f]
   *
   * @tparam Tpnu The complex type of the order
   * @tparam Tp The complex type of the argument
   * @param nu The complex order
   * @param x The complex argument
   */
  template<typename Tpnu, typename Tp>
    inline std::complex<emsr::fp_promote_t<Tpnu, Tp>>
    cyl_hankel_2(std::complex<Tpnu> nu, std::complex<Tp> x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::cyl_hankel_2<type>(nu, x);
    }

  // Spherical Hankel functions of the first kind.

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @c n
   * and complex argument @c x.
   *
   * The spherical Hankel function of the first kind is defined by
   * @f[
   *    h^{(1)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(1)}_{n+1/2}(x)
   *                 = j_n(x) + i n_n(x)
   * @f]
   * where @f$ j_n(x) @c and n_n(x) @f$ are the spherical Bessel
   * and Neumann functions respectively.
   *
   * @param n The integral order >= 0
   * @param x The complex argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    sph_hankel_1(unsigned int n, std::complex<Tp> x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_hankel_1<type>(n, x);
    }

  // Spherical Hankel functions of the second kind.

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of nonnegative order @c n
   * and complex argument @c x.
   *
   * The spherical Hankel function of the second kind is defined by
   * @f[
   *    h^{(2)}_n(x) = \left(\frac{\pi}{2x} \right) ^{1/2} H^{(2)}_{n+1/2}(x)
   *                 = j_n(x) - i n_n(x)
   * @f]
   * where @f$ j_n(x) @c and n_n(x) @f$ are the spherical Bessel
   * and Neumann functions respectively.
   *
   * @param n The integral order >= 0
   * @param x The complex argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    sph_hankel_2(unsigned int n, std::complex<Tp> x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sph_hankel_2<type>(n, x);
    }

} // namespace emsr

#endif // SF_HANKEL_H
