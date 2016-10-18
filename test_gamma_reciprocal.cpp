/*
  $HOME/bin/bin/g++ -std=gnu++17 -I. -o test_gamma_reciprocal test_gamma_reciprocal.cpp -lquadmath
  ./test_gamma_reciprocal > test_gamma_reciprocal.txt
*/

/**
 * Look at the formula for the reciprocal of the gamma for the Temme gamma
 * \frac{1}{\Gamma(1 +- \mu)}
 */

#include <bits/specfun.h>
#include <bits/float128.h>
#include <ext/math_const.h>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

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

namespace std
{
namespace __detail
{

  /**
   * 
   */
  template<typename _Tp>
    std::vector<__num_traits_t<_Tp>>
    __gamma_reciprocal_coef(std::size_t __n, _Tp __proto = _Tp{})
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
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
		   * (_Real{1} + __riemann_zeta_m_1(_Real(__j + 1 - __k)));
	  __c.push_back((_S_gamma_e * __c[__j] + __sign(__j) * __sum) / __j);
	}
      return __c;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __gamma_reciprocal(_Tp __a)
    {
      std::array<long double, 26>
      _S_c
      {{
	 0.000000000000000000e+00L,
	 1.000000000000000000e+00L,
	 5.772156649015328606e-01L,
	-6.558780715202538810e-01L,
	-4.200263503409523551e-02L,
	 1.665386113822914895e-01L,
	-4.219773455554433676e-02L,
	-9.621971527876973561e-03L,
	 7.218943246663099543e-03L,
	-1.165167591859065112e-03L,
	-2.152416741149509656e-04L,
	 1.280502823881161902e-04L,
	-2.013485478078824803e-05L,
	-1.250493482142674951e-06L,
	 1.133027231981700022e-06L,
	-2.056338416977613597e-07L,
	 6.116095104478334562e-09L,
	 5.002007644476128389e-09L,
	-1.181274570484747060e-09L,
	 1.043426711616684967e-10L,
	 7.782263434328190273e-12L,
	-3.696805613340574099e-12L,
	 5.100370285646940727e-13L,
	-2.058326339174801492e-14L,
	-5.348122461614610358e-15L,
	 1.226777080156782648e-15L,
      }};
      auto __ak = _Tp{1};
      auto __gam = _Tp{0};
      for (auto __k = 1u; __k < _S_c.size(); ++__k)
	__gam += _S_c[__k] * (__ak *= __a);
      return __gam;
    }

} // namespace __detail
} // namespace std

template<typename _Tp>
  test_gamma_reciprocal(_Tp __proto)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::size_t __n = 50;
    auto __c = std::__detail::__gamma_reciprocal_coef<_Real>(__n, __proto);

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
	      << ' ' << std::setw(width) << "1/Gamma(a)"
	      << ' ' << std::setw(width) << "1/std::tgamma(a)"
	      << ' ' << std::setw(width) << "delta"
    	      << '\n';
    for (auto __k = 0; __k < 100; ++__k)
      {
	auto __a = __k * 0.01;
	auto __gammarg = std::__detail::__gamma_reciprocal(__a);
	auto __gammars = 1 / std::tgamma(__a);
	std::cout << ' ' << std::setw(4) << __a
	        << ' ' << std::setw(width) << __gammarg
	        << ' ' << std::setw(width) << __gammars
		<< ' ' << std::setw(width) << (__gammarg - __gammars) / __gammars
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
}
