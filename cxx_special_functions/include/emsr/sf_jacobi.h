#ifndef SF_JACOBI_H
#define SF_JACOBI_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_jacobi.tcc>

namespace emsr
{

  // Jacobi polynomials

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @c n and @c float orders @f$ \alpha, \beta > -1 @f$
   * and argument @c x.
   *
   * @see jacobi for details.
   */
  inline float
  jacobif(unsigned n, float alpha, float beta, float x)
  {
    return emsr::detail::jacobi_recur<float>(n, alpha, beta, x).P_n;
  }

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @c n and @c <tt>long double</tt> orders
   * @f$ \alpha, \beta > -1 @f$ and argument @c x.
   *
   * @see jacobi for details.
   */
  inline long double
  jacobil(unsigned n, long double alpha, long double beta, long double x)
  {
    return emsr::detail::jacobi_recur<long double>(n, alpha, beta, x).P_n;
  }

  /**
   * Return the Jacobi polynomial @f$ P_n^{(\alpha,\beta)}(x) @f$
   * of degree @c n and @c float orders @f$ \alpha, \beta > -1 @f$
   * and argument @c x.
   *
   * The Jacobi polynomials are generated by a three-term recursion relation:
   * @f[
   *   2 n(\alpha + \beta + n) (\alpha + \beta + 2n - 2)
   *         P^{(\alpha, \beta)}_{n}(x)
   *     = (\alpha + \beta + 2n - 1)
   *       [(\alpha^2 - \beta^2)
   *        + x(\alpha + \beta + 2n - 2)(\alpha + \beta + 2n)]
   *         P^{(\alpha, \beta)}_{n-1}(x)
   *     - 2 (\alpha + n - 1)(\beta + n - 1)(\alpha + \beta + 2n)
   *         P^{(\alpha, \beta)}_{n-2}(x)
   * @f]
   * where @f$ P_0^{(\alpha,\beta)}(x) = 1 @f$ and
   * @f$ P_1^{(\alpha,\beta)}(x)
   *      = [(\alpha - \beta) + (\alpha + \beta + 2) x] / 2 @f$.
   *
   * @tparam Talpha The real type of the order @f$ \alpha @f$
   * @tparam Tbeta The real type of the order @f$ \beta @f$
   * @tparam Tp The real type of the argument
   * @param n The non-negative integral degree
   * @param alpha The real order
   * @param beta The real order
   * @param x The real argument
   */
  template<typename Talpha, typename Tbeta, typename Tp>
    inline emsr::fp_promote_t<Talpha, Tbeta, Tp>
    jacobi(unsigned n, Talpha alpha, Tbeta beta, Tp x)
    {
      using type = emsr::fp_promote_t<Talpha, Tbeta, Tp>;
      return emsr::detail::jacobi_recur<type>(n, alpha, beta, x).P_n;
    }

  // Zernike polynomials

  /**
   * Return the Zernike polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @c n, signed order @c m,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * @see zernike for details.
   */
  inline float
  zernikef(unsigned int n, int m, float rho, float phi)
  { return emsr::detail::zernike<float>(n, m, rho, phi); }

  /**
   * Return the Zernike polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @c n, signed order @c m,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * @see zernike for details.
   */
  inline long double
  zernikel(unsigned int n, int m, long double rho, long double phi)
  { return emsr::detail::zernike<long double>(n, m, rho, phi); }

  /**
   * Return the Zernike polynomial @f$ Z_n^m(\rho,\phi) @f$
   * for non-negative degree @c n, signed order @c m,
   * and real radial argument @f$ \rho @f$ and azimuthal angle @f$ \phi @f$.
   *
   * The even Zernike polynomials are defined by:
   * @f[
   *    Z_n^m(\rho,\phi) = R_n^m(\rho)\cos(m\phi)
   * @f]
   * and the odd Zernike polynomials are defined by:
   * @f[
   *    Z_n^{-m}(\rho,\phi) = R_n^m(\rho)\sin(m\phi)
   * @f]
   * for non-negative degree @c m and @f$ m <= n @f$
   * and where @f$ R_n^m(\rho) @f$ is the radial polynomial (see radpoly).
   *
   * @tparam _Trho The real type of the radial coordinate
   * @tparam _Tphi The real type of the azimuthal angle
   * @param n The non-negative degree.
   * @param m The (signed) azimuthal order
   * @param rho The radial coordinate
   * @param phi The azimuthal angle
   */
  template<typename _Trho, typename _Tphi>
    inline emsr::fp_promote_t<_Trho, _Tphi>
    zernike(unsigned int n, int m, _Trho rho, _Tphi phi)
    {
      using type = emsr::fp_promote_t<_Trho, _Tphi>;
      return emsr::detail::zernike<type>(n, m, rho, phi);
    }

  // Radial polynomials

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @c n, order @f$ m <= n @f$, and @c float radial
   * argument @f$ \rho @f$.
   *
   * @see radpoly for details.
   */
  inline float
  radpolyf(unsigned int n, unsigned int m, float rho)
  { return emsr::detail::radial_jacobi(n, m, rho); }

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @c n, order @f$ m <= n @f$, and <tt>long double</tt> radial
   * argument @f$ \rho @f$.
   *
   * @see radpoly for details.
   */
  inline long double
  radpolyl(unsigned int n, unsigned int m, long double rho)
  { return emsr::detail::radial_jacobi(n, m, rho); }

  /**
   * Return the radial polynomial @f$ R_n^m(\rho) @f$ for non-negative
   * degree @c n, order @f$ m <= n @f$, and real radial
   * argument @f$ \rho @f$.
   *
   * The radial polynomials are defined by 
   * @f[
   *     R_n^m(\rho) = \sum_{k=0}^{\frac{n-m}{2}}
   *       \frac{(-1)^k(n-k)!}{k!(\frac{n+m}{2}-k)!(\frac{n-m}{2}-k)!}
   *       \rho^{n-2k}
   * @f]
   * for @f$ n - m @f$ even and identically 0 for @f$ n - m @f$ odd.
   * The radial polynomials can be related to the Jacobi polynomials:
   * @f[
   *    R_n^m(\rho) = (-1)^{(n-m)/2} \rho^m P_{(n-m)/2}^{(m,0)}
   * @f]
   * @see jacobi for details on the Jacobi polynomials (see jacobi).
   *
   * @tparam Tp The real type of the radial coordinate
   * @param n The non-negative degree.
   * @param m The non-negative azimuthal order
   * @param rho The radial argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    radpoly(unsigned int n, unsigned int m, Tp rho)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::radial_jacobi<type>(n, m, rho);
    }

} // namespace emsr

#endif // SF_JACOBI_H