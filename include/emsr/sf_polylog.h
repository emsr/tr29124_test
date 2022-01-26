#ifndef SF_POLYLOG_H
#define SF_POLYLOG_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/sf_polylog.tcc>

namespace emsr
{

  // Polylogarithm functions

  /**
   * Return the real polylogarithm function of real order @c s
   * and real argument @c w.
   *
   * @see polylog for details.
   */
  inline float
  polylogf(float s, float w)
  { return emsr::detail::polylog<float>(s, w); }

  /**
   * Return the complex polylogarithm function of real order @c s
   * and argument @c w.
   *
   * @see polylog for details.
   */
  inline long double
  polylogl(long double s, long double w)
  { return emsr::detail::polylog<long double>(s, w); }

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
  template<typename _Tp, typename _Wp>
    inline emsr::fp_promote_t<_Tp, _Wp>
    polylog(_Tp s, _Wp w)
    {
      using type = emsr::fp_promote_t<_Tp, _Wp>;
      return emsr::detail::polylog<type>(s, w);
    }

  /**
   * Return the complex polylogarithm function of real order @c s
   * and complex argument @c w.
   *
   * @see polylog for details.
   */
  inline std::complex<float>
  polylogf(float s, std::complex<float> w)
  { return emsr::detail::polylog<float>(s, w); }

  /**
   * Return the complex polylogarithm function of real order @c s
   * and complex argument @c w.
   *
   * @see polylog for details.
   */
  inline std::complex<long double>
  polylogl(long double s, std::complex<long double> w)
  { return emsr::detail::polylog<long double>(s, w); }

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
  template<typename _Tp, typename _Wp>
    inline std::complex<emsr::fp_promote_t<_Tp, _Wp>>
    polylog(_Tp s, std::complex<_Tp> w)
    {
      using type = emsr::fp_promote_t<_Tp, _Wp>;
      return emsr::detail::polylog<type>(s, w);
    }

  // Dirichlet eta function

  /**
   * Return the Dirichlet eta function of real argument @c s.
   *
   * @see dirichlet_eta for details.
   */
  inline float
  dirichlet_etaf(float s)
  { return emsr::detail::dirichlet_eta<float>(s); }

  /**
   * Return the Dirichlet eta function of real argument @c s.
   *
   * @see dirichlet_eta for details.
   */
  inline long double
  dirichlet_etal(long double s)
  { return emsr::detail::dirichlet_eta<long double>(s); }

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
  template<typename _Tp>
    inline _Tp
    dirichlet_eta(_Tp s)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::dirichlet_eta<type>(s);
    }

  // Dirichlet beta function

  /**
   * Return the Dirichlet beta function of real argument @c s.
   *
   * @see dirichlet_beta for details.
   */
  inline float
  dirichlet_betaf(float s)
  { return emsr::detail::dirichlet_beta<float>(s); }

  /**
   * Return the Dirichlet beta function of real argument @c s.
   *
   * @see dirichlet_beta for details.
   */
  inline long double
  dirichlet_betal(long double s)
  { return emsr::detail::dirichlet_beta<long double>(s); }

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
  template<typename _Tp>
    inline _Tp
    dirichlet_beta(_Tp s)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::dirichlet_beta<type>(s);
    }

  // Dirichlet lambda function

  /**
   * Return the Dirichlet lambda function of real argument @c s.
   *
   * @see dirichlet_lambda for details.
   */
  inline float
  dirichlet_lambdaf(float s)
  { return emsr::detail::dirichlet_lambda<float>(s); }

  /**
   * Return the Dirichlet lambda function of real argument @c s.
   *
   * @see dirichlet_lambda for details.
   */
  inline long double
  dirichlet_lambdal(long double s)
  { return emsr::detail::dirichlet_lambda<long double>(s); }

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
  template<typename _Tp>
    inline _Tp
    dirichlet_lambda(_Tp s)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::dirichlet_lambda<type>(s);
    }

  // Clausen Sl functions

  /**
   * Return the Clausen sine function @f$ Sl_m(x) @f$ of order @c m
   * and @c float argument @c x.
   *
   * @see clausen_sl for details.
   */
  inline float
  clausen_slf(unsigned int m, float x)
  { return emsr::detail::clausen_sl<float>(m, x); }

  /**
   * Return the Clausen sine function @f$ Sl_m(x) @f$ of order @c m
   * and <tt>long double</tt> argument @c x.
   *
   * @see clausen_sl for details.
   */
  inline long double
  clausen_sll(unsigned int m, long double x)
  { return emsr::detail::clausen_sl<long double>(m, x); }

  /**
   * Return the Clausen sine function @f$ Sl_m(x) @f$ of order @c m
   * and real argument @c x.
   *
   * The Clausen sine function is defined by
   * @f[
   *    Sl_m(x) = \sum_{k=1}^\infty\frac{\sin(kx)}{k^m}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param m The unsigned integer order
   * @param x The real argument
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    clausen_sl(unsigned int m, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::clausen_sl<type>(m, x);
    }

  // Clausen Cl functions

  /**
   * Return the Clausen cosine function @f$ Cl_m(x) @f$ of order @c m
   * and @c float argument @c x.
   *
   * @see clausen_cl for details.
   */
  inline float
  clausen_clf(unsigned int m, float x)
  { return emsr::detail::clausen_cl<float>(m, x); }

  /**
   * Return the Clausen cosine function @f$ Cl_m(x) @f$ of order @c m
   * and <tt>long double</tt> argument @c x.
   *
   * @see clausen_cl for details.
   */
  inline long double
  clausen_cll(unsigned int m, long double x)
  { return emsr::detail::clausen_cl<long double>(m, x); }

  /**
   * Return the Clausen cosine function @f$ Cl_m(x) @f$ of order @c m
   * and real argument @c x.
   *
   * The Clausen cosine function is defined by
   * @f[
   *    Cl_m(x) = \sum_{k=1}^\infty\frac{\cos(kx)}{k^m}
   * @f]
   *
   * @tparam _Tp The real type of the argument
   * @param m The unsigned integer order
   * @param x The real argument
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    clausen_cl(unsigned int m, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::clausen_cl<type>(m, x);
    }

  // Clausen functions - real argument

  /**
   * Return the Clausen function @f$ C_m(x) @f$ of integer order @c m
   * and @c float argument @c x.
   *
   * @see clausen for details.
   */
  inline float
  clausenf(unsigned int m, float x)
  { return emsr::detail::clausen<float>(m, x); }

  /**
   * Return the Clausen function @f$ C_m(x) @f$ of integer order @c m
   * and <tt>long double</tt> argument @c x.
   *
   * @see clausen for details.
   */
  inline long double
  clausenl(unsigned int m, long double x)
  { return emsr::detail::clausen<long double>(m, x); }

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
   * @tparam _Tp The real type of the argument
   * @param m The integral order
   * @param x The real argument
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    clausen(unsigned int m, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::clausen<type>(m, x);
    }

  // Clausen functions - complex argument

  /**
   * Return the Clausen function @f$ C_m(z) @f$ of integer order @c m
   * and <tt>std::complex<float></tt> argument @c z.
   *
   * @see clausen for details.
   */
  inline std::complex<float>
  clausenf(unsigned int m, std::complex<float> z)
  { return emsr::detail::clausen<float>(m, z); }

  /**
   * Return the Clausen function @f$ C_m(z) @f$ of integer order @c m
   * and <tt>std::complex<long double></tt> argument @c z.
   *
   * @see clausen for details.
   */
  inline std::complex<long double>
  clausenl(unsigned int m, std::complex<long double> z)
  { return emsr::detail::clausen<long double>(m, z); }

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
   * @tparam _Tp The real type of the complex components
   * @param m The integral order
   * @param z The complex argument
   */
  template<typename _Tp>
    inline std::complex<emsr::fp_promote_t<_Tp>>
    clausen(unsigned int m, std::complex<_Tp> z)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::clausen<type>(m, z);
    }

  // Fermi-Dirac integrals.

  /**
   * Return the Fermi-Dirac integral of @c float order s and argument x.
   *
   * @see fermi_dirac for details.
   */
  inline float
  fermi_diracf(float s, float x)
  { return emsr::detail::fermi_dirac<float>(s, x); }

  /**
   * Return the Fermi-Dirac integral of <tt> long double </tt>
   * order s and argument x.
   *
   * @see fermi_dirac for details.
   */
  inline long double
  fermi_diracl(long double s, long double x)
  { return emsr::detail::fermi_dirac<long double>(s, x); }

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
  template<typename _Tps, typename _Tp>
    inline emsr::fp_promote_t<_Tps, _Tp>
    fermi_dirac(_Tps s, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tps, _Tp>;
      return emsr::detail::fermi_dirac<type>(s, x);
    }

  // Bose-Einstein integrals.

  /**
   * Return the Bose-Einstein integral of @c float order s and argument x.
   *
   * @see bose_einstein for details.
   */
  inline float
  bose_einsteinf(float s, float x)
  { return emsr::detail::bose_einstein<float>(s, x); }

  /**
   * Return the Bose-Einstein integral of <tt> long double </tt>
   * order s and argument x.
   *
   * @see bose_einstein for details.
   */
  inline long double
  bose_einsteinl(long double s, long double x)
  { return emsr::detail::bose_einstein<long double>(s, x); }

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
  template<typename _Tps, typename _Tp>
    inline emsr::fp_promote_t<_Tps, _Tp>
    bose_einstein(_Tps s, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tps, _Tp>;
      return emsr::detail::bose_einstein<type>(s, x);
    }

  /**
   * Return the Hurwitz zeta function of real order @c s,
   * and complex parameter @c a.
   *
   * @see hurwitz_zeta for details.
   */
  template<typename _Tp, typename _Up>
    inline std::complex<_Tp>
    hurwitz_zeta(_Tp s, std::complex<_Up> a)
    {
      using type = emsr::fp_promote_t<_Tp, _Up>;
      return emsr::detail::hurwitz_zeta_polylog<type>(s, a);
    }

  // Periodic zeta functions

  /**
   * Return the periodic zeta function of @c float argument @c x,
   * and parameter @c s.
   *
   * @see periodic_zeta for details.
   */
  inline std::complex<float>
  periodic_zetaf(float x, float s)
  { return emsr::detail::periodic_zeta<float>(x, s); }

  /**
   * Return the periodic zeta function of <tt>long double</tt>
   * argument @c x, and parameter @c s.
   *
   * @see periodic_zeta for details.
   */
  inline std::complex<long double>
  periodic_zetal(long double x, long double s)
  { return emsr::detail::periodic_zeta<long double>(x, s); }

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
  template<typename _Tp, typename _Up>
    inline emsr::fp_promote_t<std::complex<_Tp>, _Up>
    periodic_zeta(_Tp x, _Up s)
    {
      using type = emsr::fp_promote_t<_Tp, _Up>;
      return emsr::detail::periodic_zeta<type>(x, s);
    }

  /**
   * Return the periodic zeta function of complex argument @c z,
   * and real parameter @c s.
   *
   * @see periodic_zeta for details.
   */
  template<typename _Tp, typename _Up>
    inline emsr::fp_promote_t<std::complex<_Tp>, std::complex<_Up>>
    periodic_zeta(std::complex<_Up> z, _Tp s)
    {
      using type = emsr::fp_promote_t<_Tp, _Up>;
      return emsr::detail::periodic_zeta<type>(z, s);
    }

} // namespace emsr

#endif // SF_POLYLOG_H
