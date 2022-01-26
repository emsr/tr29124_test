#ifndef SF_HANKEL_H
#define SF_HANKEL_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/sf_hankel.tcc>

namespace emsr
{

  // Cylindrical Hankel functions of the first kind.

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>std::complex<float></tt> order @f$ \nu @f$
   * and argument @c x.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<float>
  cyl_hankel_1f(std::complex<float> nu, std::complex<float> x)
  { return emsr::detail::cyl_hankel_1<float>(nu, x); }

  /**
   * Return the complex cylindrical Hankel function of the first kind
   * @f$ H^{(1)}_\nu(x) @f$ of <tt>std::complex<long double></tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see cyl_hankel_1 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_1l(std::complex<long double> nu, std::complex<long double> x)
  { return emsr::detail::cyl_hankel_1<long double>(nu, x); }

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
   * @tparam _Tpnu The complex type of the order
   * @tparam _Tp The complex type of the argument
   * @param nu The complex order
   * @param x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tpnu, _Tp>>
    cyl_hankel_1(std::complex<_Tpnu> nu, std::complex<_Tp> x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::cyl_hankel_1<type>(nu, x);
    }

  // Cylindrical Hankel functions of the second kind.

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>std::complex<float></tt> order @f$ \nu @f$
   * and argument @c x.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<float>
  cyl_hankel_2f(std::complex<float> nu, std::complex<float> x)
  { return emsr::detail::cyl_hankel_2<float>(nu, x); }

  /**
   * Return the complex cylindrical Hankel function of the second kind
   * @f$ H^{(2)}_\nu(x) @f$ of <tt>std::complex<long double></tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see cyl_hankel_2 for more details.
   */
  inline std::complex<long double>
  cyl_hankel_2l(std::complex<long double> nu, std::complex<long double> x)
  { return emsr::detail::cyl_hankel_2<long double>(nu, x); }

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
   * @tparam _Tpnu The complex type of the order
   * @tparam _Tp The complex type of the argument
   * @param nu The complex order
   * @param x The complex argument
   */
  template<typename _Tpnu, typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tpnu, _Tp>>
    cyl_hankel_2(std::complex<_Tpnu> nu, std::complex<_Tp> x)
    {
      using type = emsr::fp_promote_t<_Tpnu, _Tp>;
      return emsr::detail::cyl_hankel_2<type>(nu, x);
    }

  // Spherical Hankel functions of the first kind.

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @c n
   * and <tt>std::complex<float></tt> argument @c x.
   *
   * @see sph_hankel_1 for more details.
   */
  inline std::complex<float>
  sph_hankel_1f(unsigned int n, std::complex<float> x)
  { return emsr::detail::sph_hankel_1<float>(n, x); }

  /**
   * Return the complex spherical Hankel function of the first kind
   * @f$ h^{(1)}_n(x) @f$ of non-negative integral @c n
   * and <tt>std::complex<long double></tt> argument @c x.
   *
   * @see sph_hankel_1 for more details.
   */
  inline std::complex<long double>
  sph_hankel_1l(unsigned int n, std::complex<long double> x)
  { return emsr::detail::sph_hankel_1<long double>(n, x); }

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
  template<typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tp>>
    sph_hankel_1(unsigned int n, std::complex<_Tp> x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::sph_hankel_1<type>(n, x);
    }

  // Spherical Hankel functions of the second kind.

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of non-negative integral @c n
   * and <tt>std::complex<float></tt> argument @c x.
   *
   * @see sph_hankel_2 for more details.
   */
  inline std::complex<float>
  sph_hankel_2f(unsigned int n, std::complex<float> x)
  { return emsr::detail::sph_hankel_2<float>(n, x); }

  /**
   * Return the complex spherical Hankel function of the second kind
   * @f$ h^{(2)}_n(x) @f$ of non-negative integral @c n
   * and <tt>std::complex<long double></tt> argument @c x.
   *
   * @see sph_hankel_2 for more details.
   */
  inline std::complex<long double>
  sph_hankel_2l(unsigned int n, std::complex<long double> x)
  { return emsr::detail::sph_hankel_2<long double>(n, x); }

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
  template<typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tp>>
    sph_hankel_2(unsigned int n, std::complex<_Tp> x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::sph_hankel_2<type>(n, x);
    }

} // namespace emsr

#endif // SF_HANKEL_H
