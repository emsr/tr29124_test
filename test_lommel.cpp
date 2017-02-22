/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_lommel test_lommel.cpp -lquadmath -Lwrappers -lwrap_gsl -lwrap_burkhardt
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH ./test_lommel > test_lommel.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_lommel test_lommel.cpp -lquadmath -Lwrappers -lwrap_gsl -lwrap_burkhardt
./test_lommel > test_lommel.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>

//namespace std
//{
//namespace __detail
//{

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
      if (__isnan(__mu) || __isnan(__nu) || __isnan(__z))
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
      if (__isnan(__mu) || __isnan(__nu) || __isnan(__z))
	return _S_NaN;
      else if (__nu < _Tp{0})
	return __lommel_2(__mu, -__nu, __z);
      else
	{
	  const auto im = __gnu_cxx::__fp_is_odd_integer(__mu - __nu);
	  const auto ip = __gnu_cxx::__fp_is_odd_integer(__mu + __nu);
	  if (im && im() < 0)
            {
	      return 0;
	    }
	  else if (ip && ip() < 0)
            {
	      return 0;
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
		    * (__sc.__sin_value * __Bess.__J_value
		     - __sc.__cos_value * __Bess.__N_value);
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

//} // namespace __detail
//} // namespace std

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
		<< ' ' << std::setw(width) << __lommel_1(mu, nu, z)
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
		<< ' ' << std::setw(width) << __lommel_2(mu, nu, z)
		<< '\n';
    }
  }

int
main()
{
  test_lommel_1(1.0);
  test_lommel_2(1.0);
}
