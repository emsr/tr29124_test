/*
$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_gamma_reciprocal test_gamma_reciprocal.cpp -lquadmath -lmpfr
./test_gamma_reciprocal > test_gamma_reciprocal.txt
*/

/**
 * Look at the formula for the reciprocal of the gamma for the Temme gamma
 * \frac{1}{\Gamma(1 +- \mu)}
 */

#include <bits/specfun.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include <bits/float128.h>
#include <bits/complex128.h>

#include <mpreal.h>
#include <ext/math_const.h>
#include <ext/math_const_mpreal.h>

// I'm not sure why I need this here and not other places...
template<>
  constexpr std::array<float, 7>
  std::__detail::_GammaSpouge<float>::_S_cheby;
template<>
  constexpr std::array<double, 18>
  std::__detail::_GammaSpouge<double>::_S_cheby;
template<>
  constexpr std::array<long double, 22>
  std::__detail::_GammaSpouge<long double>::_S_cheby;

template<>
  constexpr std::array<float, 7>
  std::__detail::_GammaLanczos<float>::_S_cheby;
template<>
  constexpr std::array<double, 10>
  std::__detail::_GammaLanczos<double>::_S_cheby;
template<>
  constexpr std::array<long double, 11>
  std::__detail::_GammaLanczos<long double>::_S_cheby;


  /**
   * 
   */
  template<typename _Tp>
    std::vector<std::__detail::__num_traits_t<_Tp>>
    __gamma_reciprocal_series_coef(std::size_t __n, _Tp __proto = _Tp{})
    {
      using _Val = _Tp;
      using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_gamma_e = __gnu_cxx::__const_gamma_e(std::real(__proto));
      auto __sign = [](std::size_t __i){ return __i & 1u == 1u ? -1 : +1; };
      std::vector<_Real> __c;
      __c.push_back(_Real{0});
      __c.push_back(_Real{1});
      for (auto __j = 1u; __j < __n; ++__j)
	{
	  auto __sum = _Real{0};
	  for (auto __k = 1u; __k < __j; ++__k)
	    __sum += __sign(__k) * __c[__k]
		   * (_Real{1} + std::__detail::__riemann_zeta_m_1(_Real(__j + 1 - __k)));
	  __c.push_back((_S_gamma_e * __c[__j] + __sign(__j) * __sum) / __j);
	}
      return __c;
    }

  /**
   * Return the reciprocal of the Gamma function by series.
   * The reciprocal of the Gamma function is given by
   * @f[
   *   \frac{1}{\Gamma(a)} = \sum_{k=1}^{\infty} c_k a^k
   * @f]
   * where the coefficients are defined by recursion:
   * @f[
   *   c_{k+1} = \frac{1}{k}\left[\gamma_E c_k
   *           + (-1)^k\sum_{j=1}^{k-1}(-1)^j\zeta(j+1-k)c_j\right]
   * @f]
   * where @f$ c_1 = 1 @f$
   */
  template<typename _Tp>
    _Tp
    __gamma_reciprocal_series(_Tp __a)
    {
      static constexpr std::array<long double, 31>
      _S_c
      {{
	 0.0000000000000000000000000000000000000000L,
	 1.0000000000000000000000000000000000000000L,
	 0.5772156649015328606065120900824024310422L,
	-0.6558780715202538810770195151453904812798L,
	-0.0420026350340952355290039348754298187114L,
	 0.1665386113822914895017007951021052357178L,
	-0.0421977345555443367482083012891873913017L,
	-0.0096219715278769735621149216723481989754L,
	 0.0072189432466630995423950103404465727099L,
	-0.0011651675918590651121139710840183886668L,
	-0.0002152416741149509728157299630536478065L,
	 0.0001280502823881161861531986263281643234L,
	-0.0000201348547807882386556893914210218184L,
	-0.0000012504934821426706573453594738330922L,
	 0.0000011330272319816958823741296203307449L,
	-0.0000002056338416977607103450154130020573L,
	 0.0000000061160951044814158178624986828553L,
	 0.0000000050020076444692229300556650480600L,
	-0.0000000011812745704870201445881265654365L,
	 0.0000000001043426711691100510491540332312L,
	 0.0000000000077822634399050712540499373114L,
	-0.0000000000036968056186422057081878158781L,
	 0.0000000000005100370287454475979015481323L,
	-0.0000000000000205832605356650678322242954L,
	-0.0000000000000053481225394230179823700173L,
	 0.0000000000000012267786282382607901588938L,
	-0.0000000000000001181259301697458769513765L,
	 0.0000000000000000011866922547516003325798L,
	 0.0000000000000000014123806553180317815558L,
	-0.0000000000000000002298745684435370206592L,
	 0.0000000000000000000171440632192733743338L,
      }};
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__a));
      auto __ak = _Tp{1};
      auto __gam = _Tp{0};
      for (auto __k = 1u; __k < _S_c.size(); ++__k)
	{
	  __ak *= __a;
	  auto __term = _S_c[__k] * __ak;
	  __gam += __term;
	  if (std::abs(__term) < _S_eps)
	    break;
	}
      return __gam;
    }

  /**
   * Return the reciprocal of the Gamma function by infinite product.
   * The reciprocal of the Gamma function is given by
   * @f[
   *   \frac{1}{\Gamma(a)} = ae^{\gamma_E a}\Pi_{k=1}^{\infty}
   *                     (\left 1+\frac{a}{k}\right)e^{-a/k}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __gamma_reciprocal_prod(_Tp __a)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__a));
      const auto _S_gamma_e = __gnu_cxx::__const_gamma_e(std::real(__a));
      const auto _S_max_iter = 1000;
      auto __gam = __a * std::exp(_S_gamma_e * __a);
      for (auto __k = 1u; __k < _S_max_iter; ++__k)
	{
	  auto __rat = __a / __k;
	  __gam *= (_Tp{1} + __rat) * std::exp(-__a / __k);
	  if (std::abs(__rat) < _S_eps)
	    break;
	}
      return __gam;
    }

template<typename _Tp>
  void
  test_gamma_reciprocal(_Tp __proto)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(__gnu_cxx::__digits10(__proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::size_t __n = 50;
    auto __c = __gamma_reciprocal_series_coef<_Real>(__n, __proto);

    std::cout << '\n'
	      << ' ' << std::setw(4) << "k"
	      << ' ' << std::setw(width) << "c"
    	      << '\n';
    for (auto __k = 0u; __k < __c.size(); ++__k)
      std::cout << ' ' << std::setw(4) << __k
	        << ' ' << std::setw(width) << __c[__k]
    	        << '\n';

    std::cout << '\n'
	      << ' ' << std::setw(4) << "a"
	      << ' ' << std::setw(width) << "1/Gamma(a) ser"
	      << ' ' << std::setw(width) << "1/Gamma(a) prod"
	      << ' ' << std::setw(width) << "1/std::tgamma(a)"
	      << ' ' << std::setw(width) << "delta series"
	      << ' ' << std::setw(width) << "delta product"
    	      << '\n';
    for (auto __k = 0; __k < 100; ++__k)
      {
	auto __a = __k * _Tp{0.01L};
	auto __gammargs = __gamma_reciprocal_series(__a);
	auto __gammargp = __gamma_reciprocal_prod(__a);
	auto __gammars = 1 / std::tgamma(__a);
	std::cout << ' ' << std::setw(4) << __a
	        << ' ' << std::setw(width) << __gammargs
	        << ' ' << std::setw(width) << __gammargp
	        << ' ' << std::setw(width) << __gammars
		<< ' ' << std::setw(width) << (__gammargs - __gammars) / __gammars
		<< ' ' << std::setw(width) << (__gammargp - __gammars) / __gammars
    	        << '\n';
      }
  }

int
main()
{
  test_gamma_reciprocal(1.0f);

  test_gamma_reciprocal(1.0);

  test_gamma_reciprocal(1.0l);

  //test_gamma_reciprocal(1.0q);

  //test_gamma_reciprocal(mpfr::mpreal(1,128));
}
