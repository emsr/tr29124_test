#ifndef SF_LEGENDRE_H
#define SF_LEGENDRE_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/specfun_state.h>
#include <emsr/detail/sf_legendre.tcc>

namespace emsr
{

  // Associated Legendre functions

  /**
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @c l, order @c m, and @c float argument @c x.
   *
   * @see assoc_legendre for more details.
   */
  inline float
  assoc_legendref(unsigned int l, unsigned int m, float x)
  { return emsr::detail::assoc_legendre_p<float>(l, m, x).P_lm; }

  /**
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @c l, order @c m, and @c <tt>long double</tt>
   * argument @c x.
   *
   * @see assoc_legendre for more details.
   */
  inline long double
  assoc_legendrel(unsigned int l, unsigned int m, long double x)
  { return emsr::detail::assoc_legendre_p<long double>(l, m, x).P_lm; }

  /**
   * Return the associated Legendre function @f$ P_l^m(x) @f$
   * of degree @c l, order @c m, and real argument @c x.
   *
   * The associated Legendre function is derived from the Legendre function
   * @f$ P_l(x) @f$ by the Rodrigues formula:
   * @f[
   *   P_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
   * @f]
   * @see legendre for details of the Legendre function of degree @c l
   * @note @f$ P_l^m(x) = 0 @c if m > l @f$.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  l  The degree <tt>l >= 0</tt>.
   * @param  m  The order <tt>m</tt>.
   * @param  x  The argument, <tt>abs(x) <= 1</tt>.
   * @throw std::domain_error if <tt>abs(x) > 1</tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    assoc_legendre(unsigned int l, unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::assoc_legendre_p<type>(l, m, x).P_lm;
    }

  // Legendre polynomials

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @c l and @c float argument @f$ |x| <= 0 @f$.
   *
   * @see legendre for more details.
   */
  inline float
  legendref(unsigned int l, float x)
  { return emsr::detail::legendre_p<float>(l, x).P_l; }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @c l and <tt>long double</tt> argument @f$ |x| <= 0 @f$.
   *
   * @see legendre for more details.
   */
  inline long double
  legendrel(unsigned int l, long double x)
  { return emsr::detail::legendre_p<long double>(l, x).P_l; }

  /**
   * Return the Legendre polynomial @f$ P_l(x) @f$ of nonnegative
   * degree @c l and real argument @f$ |x| <= 0 @f$.
   *
   * The Legendre function of order @c l and argument @c x,
   * @f$ P_l(x) @f$, is defined by:
   * @f[
   *   P_l(x) = \frac{1}{2^l l!}\frac{d^l}{dx^l}(x^2 - 1)^{l}
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param l The degree @f$ l >= 0 @f$
   * @param x The argument @c abs(x) <= 1
   * @throw std::domain_error if @c abs(x) > 1
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    legendre(unsigned int l, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::legendre_p<type>(l, x).P_l;
    }

  // Spherical associated Legendre functions

  /**
   * Return the spherical Legendre function of nonnegative integral
   * degree @c l and order @c m and float angle @f$ \theta @f$ in radians.
   *
   * @see sph_legendre for details.
   */
  inline float
  sph_legendref(unsigned int l, unsigned int m, float theta)
  { return detail::sph_legendre<float>(l, m, theta); }

  /**
   * Return the spherical Legendre function of nonnegative integral
   * degree @c l and order @c m and <tt>long double</tt> angle @f$ \theta @f$
   * in radians.
   *
   * @see sph_legendre for details.
   */
  inline long double
  sph_legendrel(unsigned int l, unsigned int m, long double theta)
  { return detail::sph_legendre<long double>(l, m, theta); }

  /**
   * Return the spherical Legendre function of nonnegative integral
   * degree @c l and order @c m and real angle @f$ \theta @f$
   * in radians.
   *
   * The spherical Legendre function is defined by
   * @f[
   *    Y_l^m(\theta,0) = (-1)^m\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}
   *                         P_l^m(\cos\theta)
   * @f]
   * where @f$ P_l^m(x) @f$ is the associated Legendre polynomial.
   * The full (complex) spherical harmonic function includes a phase factor
   * in the azimuthal angle @f$ \phi @f$:
   * @f[
   *    Y_l^m(\theta,\phi) = Y_l^m(\theta,0) e^{im\phi}
   * @f]
   *
   * @tparam Tp The floating-point type of the angle @c theta.
   * @param l The order <tt> l >= 0 </tt>
   * @param m The degree <tt> m >= 0 </tt> and <tt> m <= l </tt>
   * @param theta The radian polar angle argument
   * @see assoc_legendre for the unnormalized associated Legendre polynomial.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sph_legendre(unsigned int l, unsigned int m, Tp theta)
    {
      using type = emsr::fp_promote_t<Tp>;
      return detail::sph_legendre<type>(l, m, theta);
    }

  // Legendre functions of the second kind

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$
   * of nonnegative degree @c l and @c float argument.
   *
   * @see legendre_q for details.
   */
  inline float
  legendre_qf(unsigned int l, float x)
  { return emsr::detail::legendre_q<float>(l, x).Q_l; }

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$
   * of nonnegative degree @c l and <tt>long double</tt> argument.
   *
   * @see legendre_q for details.
   */
  inline long double
  legendre_ql(unsigned int l, long double x)
  { return emsr::detail::legendre_q<long double>(l, x).Q_l; }

  /**
   * Return the Legendre function of the second kind @f$ Q_l(x) @f$ of
   * nonnegative degree @c l and real argument @f$ |x| <= 0 @f$.
   *
   * The Legendre function of the second kind of order @c l
   * and argument @c x, @f$ Q_l(x) @f$, is defined by:
   * @f[
   *   Q_l(x) = \frac{1}{2} \log{\frac{x+1}{x-1}} P_l(x)
   *           - \sum_{k=0}^{l-1}\frac{(l+k)!}{(l-k)!(k!)^2 s^k}
   *             \left[\psi(l+1) - \psi(k+1)\right](x-1)^k
   * @f]
   * where @f$ P_l(x) @f$ is the Legendre polynomial of degree @c l
   * and @f$ \psi(x) @f$ is the digamma or psi function which for integral
   * argument is related to the harmonic number:
   * @f$ \psi(n) = -\gamma_E + H_n @f$.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param l The degree @f$ l >= 0 @f$
   * @param x The argument @c abs(x) <= 1
   * @throw std::domain_error if @c abs(x) > 1
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    legendre_q(unsigned int l, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::legendre_q<type>(l, x).Q_l;
    }

  // Associated Legendre functions of the second kind

  /**
   * Return the associated Legendre function @f$ Q_l^m(x) @f$
   * of degree @c l, order @c m, and @c float argument @c x.
   *
   * @see assoc_legendre_q for more details.
   */
  inline float
  assoc_legendre_qf(unsigned int l, unsigned int m, float x)
  { return emsr::detail::assoc_legendre_q<float>(l, m, x).Q_lm; }

  /**
   * Return the associated Legendre function @f$ Q_l^m(x) @f$
   * of degree @c l, order @c m, and @c <tt>long double</tt>
   * argument @c x.
   *
   * @see assoc_legendre_q for more details.
   */
  inline long double
  assoc_legendre_ql(unsigned int l, unsigned int m, long double x)
  {
    return emsr::detail::assoc_legendre_q<long double>(l, m, x).Q_lm;
  }

  /**
   * Return the associated Legendre function @f$ Q_l^m(x) @f$
   * of degree @c l, order @c m, and real argument @c x.
   *
   * The associated Legendre function is derived from the Legendre function
   * @f$ Q_l(x) @f$ by the Rodrigues formula:
   * @f[
   *   Q_l^m(x) = (1 - x^2)^{m/2}\frac{d^m}{dx^m}Q_l(x)
   * @f]
   * @see legendre for details of the Legendre function of degree @c l
   * @note @f$ Q_l^m(x) != 0 @c if m > l @f$.
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  l  The degree <tt>l >= 0</tt>.
   * @param  m  The order <tt>m</tt>.
   * @param  x  The argument, <tt>abs(x) <= 1</tt>.
   * @throw std::domain_error if <tt>abs(x) > 1</tt>.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    assoc_legendre_q(unsigned int l, unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::assoc_legendre_q<type>(l, m, x).Q_lm;
    }

  // Spherical harmonic functions

  /**
   * Return the complex spherical harmonic function of degree @c l,
   * order @c m, and @c float zenith angle @f$ \theta @f$,
   * and azimuth angle @f$ \phi @f$.
   *
   * @see sph_harmonic for details.
   */
  inline std::complex<float>
  sph_harmonicf(unsigned int l, int m,
		float theta, float phi)
  { return emsr::detail::sph_harmonic<float>(l, m, theta, phi); }

  /**
   * Return the complex spherical harmonic function of degree @c l,
   * order @c m, and <tt>long double</tt> zenith angle @f$ \theta @f$,
   * and azimuth angle @f$ \phi @f$.
   *
   * @see sph_harmonic for details.
   */
  inline std::complex<long double>
  sph_harmonicl(unsigned int l, int m,
		long double theta, long double phi)
  {
    return emsr::detail::sph_harmonic<long double>(l, m, theta, phi);
  }

  /**
   * Return the complex spherical harmonic function of degree @c l,
   * order @c m, and real zenith angle @f$ \theta @f$,
   * and azimuth angle @f$ \phi @f$.
   *
   * The spherical harmonic function is defined by:
   * @f[
   *    Y_l^m(\theta,\phi) = (-1)^m\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}
   *                     P_l^{|m|}(\cos\theta) \exp^{im\phi}
   * @f]
   * @note @f$ Y_l^m(\theta,\phi) = 0 @c if |m| > l @f$.
   *
   * @param l The order
   * @param m The degree
   * @param theta The zenith angle in radians
   * @param phi The azimuth angle in radians
   */
  template<typename _Ttheta, typename _Tphi>
    inline std::complex<emsr::fp_promote_t<_Ttheta, _Tphi>>
    sph_harmonic(unsigned int l, int m, _Ttheta theta, _Tphi phi)
    {
      using type = emsr::fp_promote_t<_Ttheta, _Tphi>;
      return emsr::detail::sph_harmonic<type>(l, m, theta, phi);
    }

} // namespace emsr

#endif // SF_LEGENDRE_H
