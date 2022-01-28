#ifndef SF_STIRLING_H
#define SF_STIRLING_H 1

#include <emsr/detail/sf_stirling.tcc>

namespace emsr
{

  /**
   * Return the Stirling number of the first kind.
   *
   * The Stirling numbers of the first kind are the coefficients of
   * the Pocchammer polynomials or the rising factorials:
   * @f[
   *   (x)_n = \sum_{k=0}^{n} \genfrac{[}{]}{0pt}{0}{n}{k} x^k
   * @f]
   *
   * The recursion is
   * @f[
   *   \genfrac{[}{]}{0pt}{0}{n+1}{m} = \genfrac{[}{]}{0pt}{0}{n}{m-1}
   *                  - n \genfrac{[}{]}{0pt}{0}{n}{m}
   * @f]
   * with starting values
   * @f[
   *   \genfrac{[}{]}{0pt}{0}{0}{0\rightarrow m} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \genfrac{[}{]}{0pt}{0}{0\rightarrow n}{0} = {1, 0, 0, ..., 0}
   * @f]
   * The Stirling number of the first kind is denoted by other symbols
   * in the literature, usually @f$ S_n^{(m)} @f$.
   */
  template<typename Tp>
    inline Tp
    stirling_1(unsigned int n, unsigned int m)
    { return emsr::detail::stirling_1<Tp>(n, m); }

  /**
   * Return a vector of Stirling numbers of the first kind.
   */
  template<typename Tp>
    inline std::vector<Tp>
    stirling_1(unsigned int n)
    { return emsr::detail::stirling_1<Tp>(n); }

  /**
   * Return the Stirling number of the second kind by series expansion
   * or by recursion.
   *
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \genfrac{\{}{\}}{0pt}{0}{n}{m}
   *      = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
   */
  template<typename Tp>
    inline Tp
    stirling_2(unsigned int n, unsigned int m)
    { return emsr::detail::stirling_2<Tp>(n, m); }

  /**
   * Return a vector of Stirling numbers of the second kind.
   */
  template<typename Tp>
    inline std::vector<Tp>
    stirling_2(unsigned int n)
    { return emsr::detail::stirling_2<Tp>(n); }

  /**
   * Return the Lah number.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    inline Tp
    lah(unsigned int n, unsigned int k)
    { return emsr::detail::lah<Tp>(n, k); }

  /**
   * Return a vector of Lah numbers.
   * Lah numbers are defined by downward recurrence:
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    lah(unsigned int n)
    { return emsr::detail::lah<Tp>(n); }

  /**
   * Return a vector of the Bell numbers
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp>
    inline std::vector<Tp>
    bell(unsigned int n)
    { return emsr::detail::bell<Tp>(n); }

  /**
   * Evaluate the Bell polynomial
   * @f[
   *   B(n,x) = \sum_{k=0}^{n}S_n^{(k)}x^k
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp, typename _Up>
    inline _Up
    bell(unsigned int n, _Up x)
    { return emsr::detail::bell(n, x); }

} // namespace emsr

#endif // SF_STIRLING_H
