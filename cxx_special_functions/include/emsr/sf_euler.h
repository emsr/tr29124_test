#ifndef SF_EULER_H
#define SF_EULER_H 1

#include <emsr/detail/sf_euler.tcc>
#include <emsr/sf_gamma.h> // binomial

namespace emsr
{

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename Tp>
    inline Tp
    euler(unsigned int n)
    { return emsr::detail::euler<Tp>(n); }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers of the first kind are defined by recursion:
   * @f[
   *   \genfrac{\langle}{\rangle}{0pt}{0}{n}{m}
   *     = (n-m)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}
   *     + (m+1)\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}
   *   \mbox{ for } n > 0
   * @f]
   * Note that @f$ A(n,m) @f$ is a common older notation.
   */
  template<typename Tp>
    inline Tp
    eulerian_1(unsigned int n, unsigned int m)
    { return emsr::detail::eulerian_1<Tp>(n, m); }

  /**
   * Return a vector of Eulerian numbers of the first kind.
   */
  template<typename Tp>
    inline std::vector<Tp>
    eulerian_1(unsigned int n)
    { return emsr::detail::eulerian_1<Tp>(n); }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n}{m}\right\rangle
   *   = (2n-m-1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m-1}\right\rangle
   *   + (m+1)
   *     \left\langle\genfrac{\langle}{\rangle}{0pt}{0}{n-1}{m}\right\rangle
   *       \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline Tp
    eulerian_2(unsigned int n, unsigned int m)
    { return emsr::detail::eulerian_2<Tp>(n, m); }

} // namespace emsr

#endif // SF_EULER_H
