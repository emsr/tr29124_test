#ifndef SF_HERMITE_H
#define SF_HERMITE_H 1

#include <emsr/fp_type_util.h>
#include <emsr/sf_hermite.tcc>

namespace emsr
{

  // Hermite polynomials

  /**
   * Return the Hermite polynomial @f$ H_n(x) @f$ of nonnegative order n
   * and float argument @c x.
   *
   * @see hermite for details.
   */
  inline float
  hermitef(unsigned int n, float x)
  { return emsr::detail::hermite<float>(n, x); }

  /**
   * Return the Hermite polynomial @f$ H_n(x) @f$ of nonnegative order n
   * and <tt>long double</tt> argument @c x.
   *
   * @see hermite for details.
   */
  inline long double
  hermitel(unsigned int n, long double x)
  { return emsr::detail::hermite<long double>(n, x); }

  /**
   * Return the Hermite polynomial @f$ H_n(x) @f$ of order n
   * and @c real argument @c x.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   * @f]
   *
   * The Hermite polynomial obeys a reflection formula:
   * @f[
   *   H_n(-x) = (-1)^n H_n(x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c x.
   * @param n The order
   * @param x The argument
   */
  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    hermite(unsigned int n, _Tp x)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return emsr::detail::hermite<type>(n, x);
    }

} // namespace emsr

#endif // SF_HERMITE_H
