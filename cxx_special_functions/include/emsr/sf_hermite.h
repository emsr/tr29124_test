#ifndef SF_HERMITE_H
#define SF_HERMITE_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_hermite.tcc>

namespace emsr
{

  // Hermite polynomials

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
   * @tparam Tp The floating-point type of the argument @c x.
   * @param n The order
   * @param x The argument
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    hermite(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::hermite<type>(n, x);
    }

} // namespace emsr

#endif // SF_HERMITE_H
