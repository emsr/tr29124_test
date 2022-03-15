#ifndef SF_BERNOULLI_H
#define SF_BERNOULLI_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_bernoulli.tcc>

namespace emsr
{

  // Bernoulli numbers

  /**
   * Return the Bernoulli number of integer order @c n.
   *
   * The Bernoulli numbers are defined by
   * @f[
   *    B_{2n} = (-1)^{n+1} 2\frac{(2n)!}{(2\pi)^{2n}} \zeta(2n),
   *    B_1 = -1/2
   * @f]
   * All odd Bernoulli numbers except @f$ B_1 @f$ are zero.
   *
   * @param n The order.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    bernoulli(unsigned int n)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::bernoulli<type>(n);
    }

  /**
   * Return the Bernoulli polynomial @f$ B_n(x) @f$ of order n at argument x.
   *
   * The values at 0 and 1 are equal to the corresponding Bernoulli number:
   * @f[
   *   B_n(0) = B_n(1) = B_n
   * @f]
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   B_n'(x) = n B_{n-1}(x)
   * @f]
   *
   * The series expansion for the Bernoulli polynomials is:
   * @f[
   *   B_n(x) = \sum_{k=0}^{n} B_k \binom{n}{k} x^{n-k}
   * @f]
   *
   * A useful argument promotion is:
   * @f[
   *   B_n(x+1) - B_n(x) = n x^{n-1}
   * @f]
   */
  template<typename Tp>
    inline Tp
    bernoulli(unsigned int n, Tp x)
    { return emsr::detail::bernoulli(n, x); }

} // namespace emsr

#endif // SF_BERNOULLI_H
