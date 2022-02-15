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

  /**
   * @brief This routine returns the Probabilists Hermite polynomial
   * 	    of order n: @f$ He_n(x) @f$ by recursion on n.
   *
   * The Probabilists Hermite polynomial is defined by:
   * @f[
   *   He_n(x) = (-1)^n e^{x^2/2} \frac{d^n}{dx^n} e^{-x^2/2}
   * @f]
   * or
   * @f[
   *   He_n(x) = \frac{1}{2^{-n/2}}H_n\left(\frac{x}{\sqrt{2}}\right)
   * @f]
   * where @f$ H_n(x) @f$ is the Physicists Hermite function.
   *
   * The Probabilists Hermite polynomial has first and second derivatives:
   * @f[
   *    He'_n(x) = n He_{n-1}(x)
   * @f]
   * and
   * @f[
   *    He''_n(x) = n(n - 1) He_{n-2}(x)
   * @f]
   *
   * The Probabilists Hermite polynomial are monic and are orthogonal
   * with respect to the weight function
   * @f[
   *   w(x) = e^{x^2/2}
   * @f]
   *
   * @param n The order of the Hermite polynomial.
   * @param x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    hermite_he(unsigned int n, Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::prob_hermite_recur<type>(n, x).He_n;
    }

} // namespace emsr

#endif // SF_HERMITE_H
