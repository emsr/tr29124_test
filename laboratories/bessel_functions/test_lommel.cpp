/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_lommel test_lommel.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH ./test_lommel > test_lommel.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -I. -o test_lommel test_lommel.cpp -lquadmath -Lwrappers/debug -lwrap_gsl -lwrap_burkhardt
./test_lommel > test_lommel.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>

namespace std
{
namespace __detail
{

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __lommel_1_series(_Tp __mu, _Tp __nu, _Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const unsigned int _S_max_iter = 100000u;
      const auto __z2 = __z * __z;
      const auto __nu2 = __nu * __nu;

      auto __term = _Tp{1};
      auto _S1 = __term;
      for (auto __k = 1u; __k < _S_max_iter; ++__k)
	{
	  const auto __mu1 = __mu + _Tp(2 * __k - 1);
	  __term *= -__z2 / (__mu1 * __mu1 - __nu2);
	  _S1 += __term;
	  if (std::abs(__term) < _S_eps * std::abs(_S1))
	    break;
	}
      _S1 *= std::pow(__z, __mu + 1);
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
  template<typename _Tp>
    inline _Tp
    __lommel_1(_Tp __mu, _Tp __nu, _Tp __z)
    {
      auto _S_NaN = __gnu_cxx::__quiet_NaN(__z);
      if (std::isnan(__mu) || std::isnan(__nu) || std::isnan(__z))
	return _S_NaN;
      else if (__nu < _Tp{0})
	return __lommel_1(__mu, -__nu, __z);
      else
	return __lommel_1_series(__mu, __nu, __z);
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
  template<typename _Tp>
    _Tp
    __lommel_2(_Tp __mu, _Tp __nu, _Tp __z)
    {
      auto _S_NaN = __gnu_cxx::__quiet_NaN(__z);
      if (std::isnan(__mu) || std::isnan(__nu) || std::isnan(__z))
	return _S_NaN;
      else if (__nu < _Tp{0})
	return __lommel_2(__mu, -__nu, __z);
      else
	{
	  const auto im = __gnu_cxx::__fp_is_odd_integer(__mu - __nu);
	  const auto ip = __gnu_cxx::__fp_is_odd_integer(__mu + __nu);
	  if (im && im() < 0)
            {
	      return _Tp{0};
	    }
	  else if (ip && ip() < 0)
            {
	      return _Tp{0};
	    }
	  else
            {
	      const auto _S1 = __lommel_1(__mu, __nu, __z);
	      const auto __sc = std::__detail::__sincos_pi((__mu - __nu) / _Tp{2});
	      const auto __Bess = std::__detail::__cyl_bessel_jn(__nu, __z);
	      const auto _S2 = _S1
			     + std::pow(_Tp{2}, __mu - 1)
			      * std::__detail::__gamma((__mu + __nu + 1)/ _Tp{2})
			      * std::__detail::__gamma((__mu - __nu + 1)/ _Tp{2})
		    * (__sc.__sin_v * __Bess.__J_value
		     - __sc.__cos_v * __Bess.__N_value);
	      return _S2;
	    }
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __lommel_2_asymp(_Tp __mu, _Tp __nu, _Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const unsigned int _S_max_iter = 100000u;
      const auto __zm2 = _Tp{1} / (__z * __z);
      const auto __nu2 = __nu * __nu;

      auto __term = _Tp{1};
      auto _S2 = __term;
      for (auto __k = 1; __k < _S_max_iter; ++__k)
	{
	  const auto __mu1 = -__mu + _Tp(2 * __k - 1);
	  __term *= -(__mu1 * __mu1 - __nu2) * __zm2;
	  _S2 += __term;
	  if (std::abs(__term) < _S_eps * std::abs(_S2))
	    break;
	}
      _S2 *= std::pow(__z, __mu - 1);
    }

} // namespace __detail
} // namespace std

namespace __gnu_cxx
{

  /**
   * Return the Lommel function of the first kind @f$ s_{\mu,nu}(z) @f$
   * for @c float arguments.
   *
   * @see lommel_1 for details.
   */
  inline float
  lommel_1f(float __mu, float __nu, float __z)
  { return std::__detail::__lommel_1<float>(__mu, __nu, __z); }

  /**
   * Return the Lommel function of the first kind @f$ s_{\mu,nu}(z) @f$
   * for <tt> long double </tt> arguments.
   *
   * @see lommel_1 for details.
   */
  inline long double
  lommel_1l(float __mu, long double __nu, long double __z)
  { return std::__detail::__lommel_1<long double>(__mu, __nu, __z); }

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
  template<typename _Tmu, typename _Tnu, typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Tmu, _Tnu, _Tp>
    lommel_1(_Tmu __mu, _Tnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::fp_promote_t<_Tmu, _Tnu, _Tp>;
      return std::__detail::__lommel_1<__type>(__mu, __nu, __z);
    }

  /**
   * Return the Lommel function of the second kind @f$ S_{\mu,nu}(z) @f$
   * for @c float arguments.
   *
   * @see lommel_2 for details.
   */
  inline float
  lommel_2f(float __mu, float __nu, float __z)
  { return std::__detail::__lommel_2<float>(__mu, __nu, __z); }

  /**
   * Return the Lommel function of the second kind @f$ S_{\mu,nu}(z) @f$
   * for <tt> long double </tt> arguments.
   *
   * @see lommel_2 for details.
   */
  inline long double
  lommel_2l(float __mu, long double __nu, long double __z)
  { return std::__detail::__lommel_2<long double>(__mu, __nu, __z); }

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
  template<typename _Tmu, typename _Tnu, typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Tmu, _Tnu, _Tp>
    lommel_2(_Tmu __mu, _Tnu __nu, _Tp __z)
    {
      using __type = __gnu_cxx::fp_promote_t<_Tmu, _Tnu, _Tp>;
      return std::__detail::__lommel_2<__type>(__mu, __nu, __z);
    }

} // namespace __gnu_cxx

template<typename _Tp>
  void
  test_lommel_1(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto mu = _Tp{6}/_Tp{5};
    auto nu = _Tp{1}/_Tp{5};
    const auto del = _Tp{1} / _Tp{10};
    std::cout << "\n\n";
    for (int i = 0; i <= +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << __gnu_cxx::lommel_1(mu, nu, z)
		<< '\n';
    }
  }

template<typename _Tp>
  void
  test_lommel_2(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto mu = _Tp{6}/_Tp{5};
    auto nu = _Tp{1}/_Tp{5};
    const auto del = _Tp{1} / _Tp{10};
    std::cout << "\n\n";
    for (int i = 0; i <= +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << __gnu_cxx::lommel_2(mu, nu, z)
		<< '\n';
    }
  }

int
main()
{
  test_lommel_1(1.0);
  test_lommel_2(1.0);
}
