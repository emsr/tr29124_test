#ifndef SF_POLYLOG_H
#define SF_POLYLOG_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_polylog.tcc>

namespace emsr
{

  // Polylogarithm functions

  /**
   * Return the polylogarithm function of real order @c s
   * and real argument @c w.
   *
   * The polylogarithm function is defined by
   * @f[
   *    Li_s(w) = \sum_{k=1}^{\infty} \frac{w^k}{k^s}
   * @f]
   *
   * @param s The order.
   * @param w Argument.
   */
  template<typename Tp, typename _Wp>
    inline emsr::fp_promote_t<Tp, _Wp>
    polylog(Tp s, _Wp w)
    {
      using type = emsr::fp_promote_t<Tp, _Wp>;
      return emsr::detail::polylog<type>(s, w);
    }

  /**
   * Return the complex polylogarithm function of real order @c s
   * and complex argument @c w.
   *
   * The polylogarithm function is defined by
   * @f[
   *    Li_s(w) = \sum_{k=1}^{\infty} \frac{w^k}{k^s}
   * @f]
   *
   * @param s The order.
   * @param w Argument.
   */
  template<typename Tp, typename _Wp>
    inline std::complex<emsr::fp_promote_t<Tp, _Wp>>
    polylog(Tp s, std::complex<Tp> w)
    {
      using type = emsr::fp_promote_t<Tp, _Wp>;
      return emsr::detail::polylog<type>(s, w);
    }

  // Dirichlet eta function

  /**
   * Return the Dirichlet eta function of real argument @c s.
   *
   * The Dirichlet eta function is defined by
   * @f[
   *    \eta(s) = \sum_{k=1}^\infty \frac{(-1)^k}{k^s}
   *    = \left( 1 - 2^{1-s} \right) \zeta(s)
   * @f]
   * An important reflection formula is:
   * @f[
   *    \eta(-s) = 2 \frac{1-2^{-s-1}}{1-2^{-s}} \pi^{-s-1} 
   *              s \sin(\frac{\pi s}{2}) \Gamma(s) \eta(s+1)
   * @f]
   * The Dirichlet eta function, in terms of the polylogarithm, is
   * @f[
   *   \eta(s) = -\Re[Li_s(-1)]
   * @f]
   *
   * @param s The order.
   */
  template<typename Tp>
    inline Tp
    dirichlet_eta(Tp s)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::dirichlet_eta<type>(s);
    }

  // Dirichlet beta function

  /**
   * Return the Dirichlet beta function of real argument @c s.
   *
   * The Dirichlet beta function is defined by:
   * @f[
   *    \beta(s) = \sum_{k=0}^\infty \frac{(-1)^k}{(2k+1)^s}
   * @f]
   * An important reflection formula is:
   * @f[
   *    \beta(1-s) = \left( \frac{2}{\pi}\right)^s \sin(\frac{\pi s}{2})
   *               \Gamma(s) \beta(s)
   * @f]
   * The Dirichlet beta function, in terms of the polylogarithm, is
   * @f[
   *   \beta(s) = \Im[Li_s(i)]
   * @f]
   *
   * @param s The order.
   */
  template<typename Tp>
    inline Tp
    dirichlet_beta(Tp s)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::dirichlet_beta<type>(s);
    }

  // Dirichlet lambda function

  /**
   * Return the Dirichlet lambda function of real argument @c s.
   *
   * The Dirichlet lambda function is defined by
   * @f[
   *    \lambda(s) = \sum_{k=0}^\infty \frac{1}{(2k+1)^s}
   *    = \left( 1 - 2^{-s} \right) \zeta(s)
   * @f]
   * In terms of the Riemann zeta and the Dirichlet eta functions
   * @f[
   *   \lambda(s) = \frac{1}{2}(\zeta(s) + \eta(s))
   * @f]
   *
   * @param s The order.
   */
  template<typename Tp>
    inline Tp
    dirichlet_lambda(Tp s)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::dirichlet_lambda<type>(s);
    }

  // Clausen Sl functions

  /**
   * Return the Clausen sine function @f$ Sl_m(x) @f$ of order @c m
   * and real argument @c x.
   *
   * The Clausen sine function is defined by
   * @f[
   *    Sl_m(x) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^m}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param m The unsigned integer order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    clausen_sl(unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::clausen_sl<type>(m, x);
    }

  // Clausen Cl functions

  /**
   * Return the Clausen cosine function @f$ Cl_m(x) @f$ of order @c m
   * and real argument @c x.
   *
   * The Clausen cosine function is defined by
   * @f[
   *    Cl_m(x) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^m}
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param m The unsigned integer order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    clausen_cl(unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::clausen_cl<type>(m, x);
    }

  // Clausen functions - real argument

  /**
   * Return the Clausen function @f$ C_m(x) @f$ of integer order @c m
   * and real argument @c x.
   *
   * The Clausen function is defined by
   * @f[
   *    C_m(x)
   *      = Sl_m(x) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^m} \mbox{ for even } m
   *      = Cl_m(x) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^m} \mbox{ for odd } m
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param m The integral order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    clausen(unsigned int m, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::clausen<type>(m, x);
    }

  // Clausen functions - complex argument

  /**
   * Return the Clausen function @f$ C_m(z) @f$ of integer order @c m
   * and complex argument @c z.
   *
   * The Clausen function is defined by
   * @f[
   *    C_m(z) = Sl_m(z) = \sum_{k=1}^\infty\frac{\sin(kz)}{k^m} \mbox{ for even } m
   *           = Cl_m(z) = \sum_{k=1}^\infty\frac{\cos(kz)}{k^m} \mbox{ for odd } m
   * @f]
   *
   * @tparam Tp The real type of the complex components
   * @param m The integral order
   * @param z The complex argument
   */
  template<typename Tp>
    inline std::complex<emsr::fp_promote_t<Tp>>
    clausen(unsigned int m, std::complex<Tp> z)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::clausen<type>(m, z);
    }

  // Fermi-Dirac integrals.

  /**
   * Return the Fermi-Dirac integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    F_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} + 1}dt
   *           = -Li_{s+1}(-e^x)
   * @f]
   *
   * @param s  The order s > -1.
   * @param x  The real argument.
   * @return  The real Fermi-Dirac integral F_s(x),
   */
  template<typename Tps, typename Tp>
    inline emsr::fp_promote_t<Tps, Tp>
    fermi_dirac(Tps s, Tp x)
    {
      using type = emsr::fp_promote_t<Tps, Tp>;
      return emsr::detail::fermi_dirac<type>(s, x);
    }

  // Bose-Einstein integrals.

  /**
   * Return the Bose-Einstein integral of integer or real order s
   * and real argument x.
   * @see https://en.wikipedia.org/wiki/Clausen_function
   * @see http://dlmf.nist.gov/25.12.16
   *
   * @f[
   *    G_s(x) = \frac{1}{\Gamma(s+1)}\int_0^\infty \frac{t^s}{e^{t-x} - 1}dt
   *           = Li_{s+1}(e^x)
   * @f]
   *
   * @param s  The order s >= 0.
   * @param x  The real argument.
   * @return  The real Bose-Einstein integral G_s(x),
   */
  template<typename Tps, typename Tp>
    inline emsr::fp_promote_t<Tps, Tp>
    bose_einstein(Tps s, Tp x)
    {
      using type = emsr::fp_promote_t<Tps, Tp>;
      return emsr::detail::bose_einstein<type>(s, x);
    }

  // Hurwitz zeta function.

  /**
   * Return the Hurwitz zeta function of real order @c s,
   * and complex parameter @c a.
   *
   * @see hurwitz_zeta for details.
   */
  template<typename Tp, typename Up>
    inline std::complex<Tp>
    hurwitz_zeta(Tp s, std::complex<Up> a)
    {
      using type = emsr::fp_promote_t<Tp, Up>;
      return emsr::detail::hurwitz_zeta_polylog<type>(s, a);
    }

  // Periodic zeta functions

  /**
   * Return the periodic zeta function of real argument @c x,
   * and parameter @c s.
   *
   * The the periodic zeta function is defined by
   * @f[
   *    F(x, s) = \sum_{n=1}^{\infty}\frac{e^{i2\pi nx}}{n^s}
   * @f]
   *
   * @param x The argument.
   * @param s The order.
   */
  template<typename Tp, typename Up>
    inline emsr::fp_promote_t<std::complex<Tp>, Up>
    periodic_zeta(Tp x, Up s)
    {
      using type = emsr::fp_promote_t<Tp, Up>;
      return emsr::detail::periodic_zeta<type>(x, s);
    }

  /**
   * Return the periodic zeta function of complex argument @c z,
   * and real parameter @c s.
   *
   * @see periodic_zeta for details.
   */
  template<typename Tp, typename Up>
    inline emsr::fp_promote_t<std::complex<Tp>, std::complex<Up>>
    periodic_zeta(std::complex<Up> z, Tp s)
    {
      using type = emsr::fp_promote_t<Tp, Up>;
      return emsr::detail::periodic_zeta<type>(z, s);
    }

} // namespace emsr

#endif // SF_POLYLOG_H
