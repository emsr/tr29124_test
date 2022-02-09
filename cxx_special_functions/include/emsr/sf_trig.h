#ifndef SF_TRIG_H
#define SF_TRIG_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_trig.tcc>

namespace emsr
{

  // Reperiodized sine function.

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized sine function is defined by:
   * @f[
   * 	\sin_\pi(x) = \sin(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sin_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sin_pi<type>(x);
    }

  /**
   * Return the reperiodized sine of complex argument z:
   * @f[
   *   \mathrm{sin_\pi}(z) = \sin(\pi z)
   *     = \mathrm{sin_\pi}(x) \mathrm{cosh_\pi}(y)
   *   + i \mathrm{cos_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename Tp>
    inline std::complex<Tp>
    sin_pi(std::complex<Tp> z)
    { return emsr::detail::sin_pi(z); }

  // Reperiodized hyperbolic sine function.

  /**
   * Return the reperiodized hyperbolic sine function @f$ \sinh_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized hyperbolic sine function is defined by:
   * @f[
   * 	\sinh_\pi(x) = \sinh(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinh_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinh_pi<type>(x);
    }

  /**
   * Return the reperiodized hyperbolic sine of complex argument z:
   * @f[
   *   \mathrm{sinh_\pi}(z) = \sinh(\pi z)
   *     = \mathrm{\sinh_\pi}(x) \mathrm{cos_\pi}(y)
   *   + i \mathrm{\cosh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename Tp>
    inline std::complex<Tp>
    sinh_pi(std::complex<Tp> z)
    { return emsr::detail::sinh_pi(z); }

  // Reperiodized cosine function.

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized cosine function is defined by:
   * @f[
   * 	\cos_\pi(x) = \cos(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    cos_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::cos_pi<type>(x);
    }

  /**
   * Return the reperiodized cosine of complex argument z:
   * @f[
   *    \mathrm{cos_\pi}(z) = \cos(\pi z)
   *       = \mathrm{cos_\pi}(x) \mathrm{cosh_\pi}(y)
   *     - i \mathrm{sin_\pi}(x) \mathrm{sinh_\pi}(y)
   * @f]
   */
  template<typename Tp>
    inline std::complex<Tp>
    cos_pi(std::complex<Tp> z)
    { return emsr::detail::cos_pi(z); }

  // Reperiodized hyperbolic cosine function.

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized hyperbolic cosine function is defined by:
   * @f[
   * 	\cosh_\pi(x) = \cosh(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    cosh_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::cosh_pi<type>(x);
    }

  /**
   * Return the reperiodized hyperbolic cosine of complex argument z:
   * @f[
   *    \mathrm{cosh_\pi}(z) = \mathrm{cosh_\pi}(z)
   *       = \mathrm{cosh_\pi}(x) \mathrm{cos_\pi}(y)
   *     + i \mathrm{sinh_\pi}(x) \mathrm{sin_\pi}(y)
   * @f]
   */
  template<typename Tp>
    inline std::complex<Tp>
    cosh_pi(std::complex<Tp> z)
    { return emsr::detail::cosh_pi(z); }

  // Reperiodized tangent function.

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized tangent function is defined by:
   * @f[
   * 	\tan_\pi(x) = \tan(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    tan_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::tan_pi<type>(x);
    }

  /**
   * Return the reperiodized tangent of complex argument z:
   * @f[
   *   \mathrm{tan_\pi}(z) = \tan(\pi z)
   *     = \frac{\mathrm{tan_\pi}(x) + i \mathrm{tanh_\pi}(y)}
   *            {1 - i \mathrm{tan_\pi}(x) \mathrm{tanh_\pi}(y)}
   * @f]
   */
  template<typename Tp>
    inline std::complex<Tp>
    tan_pi(std::complex<Tp> z)
    { return emsr::detail::tan_pi(z); }

  // Reperiodized hyperbolic tangent function.

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for real argument @c x.
   *
   * The reperiodized hyperbolic tangent function is defined by:
   * @f[
   * 	\tanh_\pi(x) = \tanh(\pi x)
   * @f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    tanh_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::tanh_pi<type>(x);
    }

  /**
   * Return the reperiodized hyperbolic tangent of complex argument z:
   * @f[
   *   \mathrm{tanh_\pi}(z) = \tanh(\pi z)
   *     = \frac{\mathrm{tanh_\pi}(x) + i \mathrm{tan_\pi}(y)}
   *            {1 + i \mathrm{tanh_\pi}(x) \mathrm{tan_\pi}(y)}
   * @f]
   */
  template<typename Tp>
    std::complex<Tp>
    tanh_pi(std::complex<Tp> z)
    { return emsr::detail::tanh_pi(z); }

  // Sincos functions.

  /**
   * Return both the sine and the cosine of a @c double argument.
   *
   * @see sincos for details.
   */
  inline sincos_t<double>
  sincos(double x)
  { return emsr::detail::sincos<double>(x); }

  /**
   * Return both the sine and the cosine of a reperiodized argument.
   * @f[
   *   sincos(x) = {\sin(x), \cos(x)}
   * @f]
   */
  template<typename Tp>
    inline sincos_t<emsr::fp_promote_t<Tp>>
    sincos(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sincos<type>(x);
    }

  // Reperiodized sincos functions.

  /**
   * Return both the sine and the cosine of a reperiodized real argument.
   *
   * @f[
   *   sincos_\pi(x) = {\sin(\pi x), \cos(\pi x)}
   * @f]
   */
  template<typename Tp>
    inline sincos_t<emsr::fp_promote_t<Tp>>
    sincos_pi(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sincos_pi<type>(x);
    }

  /**
   * Reperiodized complex constructor.
   */
  template<typename Tp>
    inline std::complex<Tp>
    polar_pi(Tp rho, Tp phi_pi)
    {
      emsr::sincos_t<Tp> sc = sincos_pi(phi_pi);
      return std::complex<Tp>(rho * sc.cos_v, rho * sc.sin_v);
    }

  /**
   * Reperiodized complex constructor.
   */
  template<typename Tp>
    inline std::complex<Tp>
    polar_pi(Tp rho, const std::complex<Tp>& phi_pi)
    {
      emsr::sincos_t<Tp> sc = sincos_pi(phi_pi);
      return std::complex<Tp>(rho * sc.cos_v, rho * sc.sin_v);
    }

} // namespace emsr

#endif // SF_TRIG_H
