/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/fp_type_util.h>
#include <emsr/sf_gamma.h>
#include <emsr/sf_bessel.h>

namespace emsr
{
namespace detail
{

  /**
   * 
   */
  template<typename Tp>
    Tp
    lommel_1_series(Tp mu, Tp nu, Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
      const auto z2 = z * z;
      const auto nu2 = nu * nu;

      auto term = Tp{1};
      auto _S1 = term;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  const auto mu1 = mu + Tp(2 * k - 1);
	  term *= -z2 / (mu1 * mu1 - nu2);
	  _S1 += term;
	  if (std::abs(term) < s_eps * std::abs(_S1))
	    break;
	}
      _S1 *= std::pow(z, mu + 1);
      return _S1;
    }

  /**
   * Return the Lommel function of the first kind.
   * @f[
   *    s_{\mu,\nu}(z) = z^{\mu+1}\sum_{k=0}^{\infty}
   *             \frac{(-1)^kz^{2k}}{a_{k+1}(\mu,\nu)}
   * @f]
   * where
   * @f[
   *    a_{k+1}(\mu,\nu) = \prod_{m=1}^{k}\left[(\mu+2m-1)^2-\nu^2\right]
   * @f]
   */
  template<typename Tp>
    inline Tp
    lommel_1(Tp mu, Tp nu, Tp z)
    {
      auto s_NaN = emsr::quiet_NaN(z);
      if (std::isnan(mu) || std::isnan(nu) || std::isnan(z))
	return s_NaN;
      else if (nu < Tp{0})
	return lommel_1(mu, -nu, z);
      else
	return lommel_1_series(mu, nu, z);
    }


  /**
   * Return the Lommel function of the second kind.
   * @f[
   *   S_{\mu,\nu}(z) = s_{\mu,\nu}(z)
   *        + 2^{\mu-1}\Gamma\left(\frac{\mu+\nu+1}{2}\right)
   *                   \Gamma\left(\frac{\mu-\nu+1}{2}\right)
   *    \left[\sin\left(\frac{(\mu-\nu)\pi}{2}\right)J_\nu(z)
   *        - \cos\left(\frac{(\mu-\nu)\pi}{2}\right)N_\nu(z0  \right]
   * @f]
   */
  template<typename Tp>
    Tp
    lommel_2(Tp mu, Tp nu, Tp z)
    {
      auto s_NaN = emsr::quiet_NaN(z);
      if (std::isnan(mu) || std::isnan(nu) || std::isnan(z))
	return s_NaN;
      else if (nu < Tp{0})
	return lommel_2(mu, -nu, z);
      else
	{
	  const auto im = emsr::fp_is_odd_integer(mu - nu);
	  const auto ip = emsr::fp_is_odd_integer(mu + nu);
	  if (im && im() < 0)
            {
	      return Tp{0};
	    }
	  else if (ip && ip() < 0)
            {
	      return Tp{0};
	    }
	  else
            {
	      const auto _S1 = lommel_1(mu, nu, z);
	      const auto sc = emsr::detail::sincos_pi((mu - nu) / Tp{2});
	      const auto Bess = emsr::detail::cyl_bessel_jn(nu, z);
	      const auto _S2 = _S1
			     + std::pow(Tp{2}, mu - 1)
			      * emsr::detail::gamma((mu + nu + 1)/ Tp{2})
			      * emsr::detail::gamma((mu - nu + 1)/ Tp{2})
		    * (sc.sin_v * Bess.J_value
		     - sc.cos_v * Bess.N_value);
	      return _S2;
	    }
	}
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    lommel_2_asymp(Tp mu, Tp nu, Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
      const auto zm2 = Tp{1} / (z * z);
      const auto nu2 = nu * nu;

      auto term = Tp{1};
      auto _S2 = term;
      for (auto k = 1; k < s_max_iter; ++k)
	{
	  const auto mu1 = -mu + Tp(2 * k - 1);
	  term *= -(mu1 * mu1 - nu2) * zm2;
	  _S2 += term;
	  if (std::abs(term) < s_eps * std::abs(_S2))
	    break;
	}
      _S2 *= std::pow(z, mu - 1);
    }

} // namespace detail
} // namespace emsr

namespace emsr
{

  /**
   * Return the Lommel function of the first kind @f$ s_{\mu,nu}(z) @f$
   * for @c float arguments.
   *
   * @see lommel_1 for details.
   */
  inline float
  lommel_1f(float mu, float nu, float z)
  { return emsr::detail::lommel_1<float>(mu, nu, z); }

  /**
   * Return the Lommel function of the first kind @f$ s_{\mu,nu}(z) @f$
   * for <tt> long double </tt> arguments.
   *
   * @see lommel_1 for details.
   */
  inline long double
  lommel_1l(float mu, long double nu, long double z)
  { return emsr::detail::lommel_1<long double>(mu, nu, z); }

  /**
   * Return the Lommel function of the first kind.
   * @f[
   *    s_{\mu,\nu}(z) = z^{\mu+1}\sum_{k=0}^{\infty}
   *             \frac{(-1)^kz^{2k}}{a_{k+1}(\mu,\nu)}
   * @f]
   * where
   * @f[
   *    a_{k+1}(\mu,\nu) = \prod_{m=1}^{k}\left[(\mu+2m-1)^2-\nu^2\right]
   * @f]
   */
  template<typename _Tmu, typename _Tnu, typename Tp>
    inline emsr::fp_promote_t<_Tmu, _Tnu, Tp>
    lommel_1(_Tmu mu, _Tnu nu, Tp z)
    {
      using type = emsr::fp_promote_t<_Tmu, _Tnu, Tp>;
      return emsr::detail::lommel_1<type>(mu, nu, z);
    }

  /**
   * Return the Lommel function of the second kind @f$ S_{\mu,nu}(z) @f$
   * for @c float arguments.
   *
   * @see lommel_2 for details.
   */
  inline float
  lommel_2f(float mu, float nu, float z)
  { return emsr::detail::lommel_2<float>(mu, nu, z); }

  /**
   * Return the Lommel function of the second kind @f$ S_{\mu,nu}(z) @f$
   * for <tt> long double </tt> arguments.
   *
   * @see lommel_2 for details.
   */
  inline long double
  lommel_2l(float mu, long double nu, long double z)
  { return emsr::detail::lommel_2<long double>(mu, nu, z); }

  /**
   * Return the Lommel function of the second kind.
   * @f[
   *   S_{\mu,\nu}(z) = s_{\mu,\nu}(z)
   *        + 2^{\mu-1}\Gamma\left(\frac{\mu+\nu+1}{2}\right)
   *                   \Gamma\left(\frac{\mu-\nu+1}{2}\right)
   *    \left[\sin\left(\frac{(\mu-\nu)\pi}{2}\right)J_\nu(z)
   *        - \cos\left(\frac{(\mu-\nu)\pi}{2}\right)N_\nu(z0  \right]
   * @f]
   * where @f$ s_{\mu,\nu}(z) @f$ is the Lommel function of the first kind
   * @see lommel_1 and @f$ J_\nu(z) @f$ is the cylindrical Bessel function
   * @see cyl_bessel_j and @f$ N_\nu(z) @f$ is the cylindrical Neumann function
   * @see cyl_neumann.
   */
  template<typename _Tmu, typename _Tnu, typename Tp>
    inline emsr::fp_promote_t<_Tmu, _Tnu, Tp>
    lommel_2(_Tmu mu, _Tnu nu, Tp z)
    {
      using type = emsr::fp_promote_t<_Tmu, _Tnu, Tp>;
      return emsr::detail::lommel_2<type>(mu, nu, z);
    }

} // namespace emsr

template<typename Tp>
  void
  test_lommel_1(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto mu = Tp{6}/Tp{5};
    auto nu = Tp{1}/Tp{5};
    const auto del = Tp{1} / Tp{10};
    std::cout << "\n\n";
    for (int i = 0; i <= +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << emsr::lommel_1(mu, nu, z)
		<< '\n';
    }
  }

template<typename Tp>
  void
  test_lommel_2(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto mu = Tp{6}/Tp{5};
    auto nu = Tp{1}/Tp{5};
    const auto del = Tp{1} / Tp{10};
    std::cout << "\n\n";
    for (int i = 0; i <= +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << emsr::lommel_2(mu, nu, z)
		<< '\n';
    }
  }

int
main()
{
  test_lommel_1(1.0);
  test_lommel_2(1.0);
}
