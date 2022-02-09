#ifndef SF_EXPINT_H
#define SF_EXPINT_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_expint.tcc>

namespace emsr
{

  // Exponential integrals

  /**
   * Return the exponential integral @f$ Ei(x) @f$ for @c real argument @c x.
   *
   * The exponential integral is given by
   * \f[
   *   Ei(x) = -\int_{-x}^\infty \frac{e^t}{t} dt
   * \f]
   *
   * @tparam Tp The floating-point type of the argument @c x.
   * @param  x  The argument of the exponential integral function.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    expint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::expint<type>(x);
    }

  //  Exponential integrals

  /**
   * Return the exponential integral @f$ E_n(x) @f$ of integral
   * order @c n and real argument @c x.
   * The exponential integral is defined by:
   * @f[
   *    E_n(x) = \int_1^\infty \frac{e^{-tx}}{t^n}dt
   * @f]
   * In particular
   * @f[
   *    E_1(x) = \int_1^\infty \frac{e^{-tx}}{t}dt = -Ei(-x)
   * @f]
   *
   * @tparam Tp The real type of the argument
   * @param n The integral order
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    expint(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::expint<type>(n, x);
    }

  // Logarithmic integrals

  /**
   * Return the logarithmic integral of argument @c x.
   *
   * The logarithmic integral is defined by
   * @f[
   *    li(x) = \int_0^x \frac{dt}{ln(t)} = Ei(ln(x))
   * @f]
   * where @f$ Ei(x) @f$ is the exponential integral (see std::expint).
   *
   * @param x The real upper integration limit
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    logint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::logint<type>(x);
    }

  // Hyperbolic sine integrals

  /**
   * Return the hyperbolic sine integral @f$ Shi(x) @f$
   * of real argument @c x.
   *
   * The hyperbolic sine integral is defined by
   * @f[
   *    Shi(x) = \int_0^x \frac{\sinh(t)}{t}dt
   * @f]
   *
   * @tparam Tp The type of the real argument
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    sinhint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::sinhint<type>(x);
    }

  // Hyperbolic cosine integrals

  /**
   * Return the hyperbolic cosine integral @f$ Chi(x) @f$
   * of real argument @c x.
   *
   * The hyperbolic cosine integral is defined by
   * @f[
   *    Chi(x) = -\int_x^\infty \frac{\cosh(t)}{t}dt
   *     = \gamma_E + ln(x) + \int_0^x \frac{\cosh(t)-1}{t}dt
   * @f]
   *
   * @tparam Tp The type of the real argument
   * @param x The real argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    coshint(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::coshint<type>(x);
    }

} // namespace emsr

#endif // SF_EXPINT_H
