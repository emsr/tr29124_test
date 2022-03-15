
#ifndef SF_THETA_H
#define SF_THETA_H 1

#include <emsr/fp_type_util.h>

namespace emsr
{
  /**
   * Slots for Jacobi elliptic function tuple.
   */
  template<typename Tp>
    struct jacobi_ellint_t
    {
      /// Jacobi sine amplitude value.
      Tp sn_value;

      /// Jacobi cosine amplitude value.
      Tp cn_value;

      /// Jacobi delta amplitude value.
      Tp dn_value;

      constexpr Tp
      am() const noexcept
      { return std::asin(sn_value); }

      constexpr Tp
      ns() const noexcept
      { return Tp{1} / sn_value; }

      constexpr Tp
      nc() const noexcept
      { return Tp{1} / cn_value; }

      constexpr Tp
      nd() const noexcept
      { return Tp{1} / dn_value; }

      constexpr Tp
      sc() const noexcept
      { return sn_value / cn_value; }

      constexpr Tp
      sd() const noexcept
      { return sn_value / dn_value; }

      constexpr Tp
      cd() const noexcept
      { return cn_value / dn_value; }

      constexpr Tp
      cs() const noexcept
      { return cn_value / sn_value; }

      constexpr Tp
      ds() const noexcept
      { return dn_value / sn_value; }

      constexpr Tp
      dc() const noexcept
      { return dn_value / cn_value; }

      constexpr Tp
      sn_deriv() const noexcept
      { return cn_value * dn_value; }

      constexpr Tp
      cn_deriv() const noexcept
      { return -sn_value * dn_value; }
    };
}

#include <emsr/detail/sf_theta.tcc>

namespace emsr
{

  // Exponential theta_1 functions.

  /**
   * Return the exponential theta-1 function @f$ \theta_1(\nu,x) @f$
   * of period @f$ \nu @f$ and argument @c x.
   *
   * The exponential theta-1 function is defined by
   * @f[
   *    \theta_1(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j - 1/2)^2}{x} \right)
   * @f]
   *
   * @param nu The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    theta_1(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::theta_1<type>(nu, x);
    }

  // Exponential theta_2 functions.

  /**
   * Return the exponential theta-2 function @f$ \theta_2(\nu,x) @f$
   * of period @f$ \nu @f$ and argument @c x.
   *
   * The exponential theta-2 function is defined by
   * @f[
   *    \theta_2(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(\nu + j)^2}{x} \right)
   * @f]
   *
   * @param nu The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    theta_2(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::theta_2<type>(nu, x);
    }

  // Exponential theta_3 functions.

  /**
   * Return the exponential theta-3 function @f$ \theta_3(\nu,x) @f$
   * of period @f$ \nu @f$ and argument @c x.
   *
   * The exponential theta-3 function is defined by
   * @f[
   *    \theta_3(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    \exp\left( \frac{-(\nu+j)^2}{x} \right)
   * @f]
   *
   * @param nu The periodic (period = 1) argument
   * @param x The argument
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    theta_3(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::theta_3<type>(nu, x);
    }

  // Exponential theta_4 functions.

  /**
   * Return the exponential theta-4 function @f$ \theta_4(\nu,x) @f$
   * of period @f$ \nu @f$ and argument @c x.
   *
   * The exponential theta-4 function is defined by
   * @f[
   *    \theta_4(\nu,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *                 \exp\left( \frac{-(\nu + j + 1/2)^2}{x} \right)
   * @f]
   *
   * @param nu The periodic (period = 1) argument
   * @param x The argument
   */
  template<typename Tpnu, typename Tp>
    inline emsr::fp_promote_t<Tpnu, Tp>
    theta_4(Tpnu nu, Tp x)
    {
      using type = emsr::fp_promote_t<Tpnu, Tp>;
      return emsr::detail::theta_4<type>(nu, x);
    }

  // Elliptic nome function

  /**
   * Return the elliptic nome function @f$ q(k) @f$ of modulus @c k.
   *
   * The elliptic nome function is defined by
   * @f[
   *    q(k) = \exp \left(-\pi\frac{K(\sqrt{1-k^2})}{K(k)} \right)
   * @f]
   * where @f$ K(k) @f$ is the complete elliptic function of the first kind.
   *
   * @tparam Tp The real type of the modulus
   * @param k The modulus @f$ -1 <= k <= +1 @f$
   */
  template<typename Tp>
    inline Tp
    ellnome(Tp k)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::ellnome<type>(k);
    }

  // Neville theta_s functions.

  /**
   * Return the Neville theta-s function @f$ \theta_s(k,x) @f$
   * of modulus @c k and argument @c x.
   *
   * The Neville theta-s function is defined by
   * @f[
   *  \theta_s(k,x) = \sqrt{\frac{\pi}{2 k k' K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   * where @f$ q(k) @f$ is the elliptic nome, @f$ K(k) @f$ is
   * the complete Legendre elliptic integral of the first kind,
   * and @f$ \theta_1(\nu,x) @f$ is the exponential theta-1 function.
   * @see ellnome, comp_ellint_1, and theta_1 for details.
   *
   * @param k The modulus @f$ -1 <= k <= +1 @f$
   * @param x The argument
   */
  template<typename Tpk, typename Tp>
    inline emsr::fp_promote_t<Tpk, Tp>
    theta_s(Tpk k, Tp x)
    {
      using type = emsr::fp_promote_t<Tpk, Tp>;
      return emsr::detail::theta_s<type>(k, x);
    }

  // Neville theta_c functions.

  /**
   * Return the Neville theta-c function @f$ \theta_c(k,x) @f$
   * of modulus @c k and argument @c x.
   *
   * The Neville theta-c function is defined by
   * @f[
   *    \theta_c(k,x) = \sqrt{\frac{\pi}{2 k K(k)}}
   *                  \theta_1\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   * where @f$ q(k) @f$ is the elliptic nome, @f$ K(k) @f$ is
   * the complete Legendre elliptic integral of the first kind,
   * and @f$ \theta_1(\nu,x) @f$ is the exponential theta-1 function.
   * @see ellnome, comp_ellint_1, and theta_1 for details.
   *
   * @param k The modulus @f$ -1 <= k <= +1 @f$
   * @param x The argument
   */
  template<typename Tpk, typename Tp>
    inline emsr::fp_promote_t<Tpk, Tp>
    theta_c(Tpk k, Tp x)
    {
      using type = emsr::fp_promote_t<Tpk, Tp>;
      return emsr::detail::theta_c<type>(k, x);
    }

  // Neville theta_d functions.

  /**
   * Return the Neville theta-d function @f$ \theta_d(k,x) @f$
   * of modulus @c k and argument @c x.
   *
   * The Neville theta-d function is defined by
   * @f[
   *    \theta_d(k,x) = \sqrt{\frac{\pi}{2K(k)}}
   *                  \theta_3\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   * where @f$ q(k) @f$ is the elliptic nome, @f$ K(k) @f$ is
   * the complete Legendre elliptic integral of the first kind,
   * and @f$ \theta_3(\nu,x) @f$ is the exponential theta-3 function.
   * @see ellnome, comp_ellint_1, and theta_3 for details.
   *
   * @param k The modulus @f$ -1 <= k <= +1 @f$
   * @param x The argument
   */
  template<typename Tpk, typename Tp>
    inline emsr::fp_promote_t<Tpk, Tp>
    theta_d(Tpk k, Tp x)
    {
      using type = emsr::fp_promote_t<Tpk, Tp>;
      return emsr::detail::theta_d<type>(k, x);
    }

  // Neville theta_n functions.

  /**
   * Return the Neville theta-n function @f$ \theta_n(k,x) @f$
   * of modulus @c k and argument @c x.
   *
   * The Neville theta-n function is defined by
   * @f[
   *  \theta_n(k,x) = \sqrt{\frac{\pi}{2k'K(k)}}
   *                  \theta_4\left(q(k),\frac{\pi x}{2K(k)}\right)
   * @f]
   * where @f$ q(k) @f$ is the elliptic nome, @f$ K(k) @f$ is
   * the complete Legendre elliptic integral of the first kind,
   * and @f$ \theta_4(\nu,x) @f$ is the exponential theta-4 function.
   * @see ellnome, comp_ellint_1, and theta_4 for details.
   *
   * @param k The modulus @f$ -1 <= k <= +1 @f$
   * @param x The argument
   */
  template<typename Tpk, typename Tp>
    inline emsr::fp_promote_t<Tpk, Tp>
    theta_n(Tpk k, Tp x)
    {
      using type = emsr::fp_promote_t<Tpk, Tp>;
      return emsr::detail::theta_n<type>(k, x);
    }

  // Jacobi theta_1 functions.

  /**
   * Return the Jacobi theta-1 function @f$ \theta_1(q,x) @f$
   * of nome @c q and argument @c x.
   *
   * The Jacobi theta-1 function is defined by
   * @f[
   *    \theta_1(q,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(q + j - 1/2)^2}{x} \right)
   * @f]
   *
   * @param q The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tpq, typename Tp>
    inline emsr::fp_promote_t<Tpq, Tp>
    jacobi_theta_1(Tpq q, Tp x)
    {
      using type = emsr::fp_promote_t<Tpq, Tp>;
      return emsr::detail::jacobi_theta_1<type>(q, x);
    }

  // Jacobi theta_2 functions.

  /**
   * Return the Jacobi theta-2 function @f$ \theta_2(q,x) @f$
   * of nome @c q and argument @c x.
   *
   * The Jacobi theta-2 function is defined by
   * @f[
   *    \theta_2(q,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    (-1)^j \exp\left( \frac{-(q + j)^2}{x} \right)
   * @f]
   *
   * @param q The periodic (period = 2) argument
   * @param x The argument
   */
  template<typename Tpq, typename Tp>
    inline emsr::fp_promote_t<Tpq, Tp>
    jacobi_theta_2(Tpq q, Tp x)
    {
      using type = emsr::fp_promote_t<Tpq, Tp>;
      return emsr::detail::jacobi_theta_2<type>(q, x);
    }

  // Jacobi theta_3 functions.

  /**
   * Return the Jacobi theta-3 function @f$ \theta_3(q,x) @f$
   * of nome @c q and argument @c x.
   *
   * The Jacobi theta-3 function is defined by
   * @f[
   *    \theta_3(q,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *    \exp\left( \frac{-(q+j)^2}{x} \right)
   * @f]
   *
   * @param q The elliptic nome
   * @param x The argument
   */
  template<typename Tpq, typename Tp>
    inline emsr::fp_promote_t<Tpq, Tp>
    jacobi_theta_3(Tpq q, Tp x)
    {
      using type = emsr::fp_promote_t<Tpq, Tp>;
      return emsr::detail::jacobi_theta_3<type>(q, x);
    }

  // Jacobi theta_4 functions.

  /**
   * Return the Jacobi theta-4 function @f$ \theta_4(q,x) @f$
   * of nome @c q and argument @c x.
   *
   * The Jacobi theta-4 function is defined by
   * @f[
   *    \theta_4(q,x) = \frac{1}{\sqrt{\pi x}} \sum_{j=-\infty}^{+\infty}
   *                 \exp\left( \frac{-(q + j + 1/2)^2}{x} \right)
   * @f]
   *
   * @param q The elliptic nome
   * @param x The argument
   */
  template<typename Tpq, typename Tp>
    inline emsr::fp_promote_t<Tpq, Tp>
    jacobi_theta_4(Tpq q, Tp x)
    {
      using type = emsr::fp_promote_t<Tpq, Tp>;
      return emsr::detail::jacobi_theta_4<type>(q, x);
    }

  // Jacobi elliptic sine amplitude functions.

  /**
   * Return the Jacobi elliptic sine amplitude function @f$ sn(k,u) @f$
   * of real modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c sn integral is defined by
   * @f[
   *    \sin(\phi) = sn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the Legendre elliptic integral
   * of the first kind (see ellint_1).
   *
   * @tparam Kp The type of the real modulus
   * @tparam Up The type of the real argument
   * @param k The real modulus
   * @param u The real argument
   */
  template<typename Kp, typename Up>
    inline emsr::fp_promote_t<Kp, Up>
    jacobi_sn(Kp k, Up u)
    {
      using type = emsr::fp_promote_t<Kp, Up>;
      return emsr::detail::jacobi_ellint<type>(k, u).sn_value;
    }

  // Jacobi elliptic cosine amplitude functions.

  /**
   * Return the Jacobi elliptic cosine amplitude function @f$ cn(k,u) @f$
   * of real modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c cn integral is defined by
   * @f[
   *    \cos(\phi) = cn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the Legendre elliptic integral
   * of the first kind (see ellint_1).
   *
   * @tparam Kp The type of the real modulus
   * @tparam Up The type of the real argument
   * @param k The real modulus
   * @param u The real argument
   */
  template<typename Kp, typename Up>
    inline emsr::fp_promote_t<Kp, Up>
    jacobi_cn(Kp k, Up u)
    {
      using type = emsr::fp_promote_t<Kp, Up>;
      return emsr::detail::jacobi_ellint<type>(k, u).cn_value;
    }

  // Jacobi elliptic delta amplitude functions.

  /**
   * Return the Jacobi elliptic delta amplitude function @f$ dn(k,u) @f$
   * of real modulus @c k and argument @c u.
   *
   * The Jacobi elliptic @c dn integral is defined by
   * @f[
   *    \sqrt{1 - k^2\sin(\phi)} = dn(k, F(k,\phi))
   * @f]
   * where @f$ F(k,\phi) @f$ is the Legendre elliptic integral
   * of the first kind (see ellint_1).
   *
   * @tparam Kp The type of the real modulus
   * @tparam Up The type of the real argument
   * @param k The real modulus
   * @param u The real argument
   */
  template<typename Kp, typename Up>
    inline emsr::fp_promote_t<Kp, Up>
    jacobi_dn(Kp k, Up u)
    {
      using type = emsr::fp_promote_t<Kp, Up>;
      return emsr::detail::jacobi_ellint<type>(k, u).dn_value;
    }

} // namespace emsr

#endif // SF_THETA_H
