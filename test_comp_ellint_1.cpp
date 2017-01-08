/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_comp_ellint_1 test_comp_ellint_1.cpp -lquadmath
./test_comp_ellint_1 > test_comp_ellint_1.txt

$HOME/bin/bin/g++ -std=gnu++17  -g -Wall -Wextra-I. -o test_comp_ellint_1 test_comp_ellint_1.cpp -lquadmath
./test_comp_ellint_1 > test_comp_ellint_1.txt

*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
#include <bits/float128_io.h>

//  Use AGM to do an ab initio calculation of K(k).
template<typename _Tp>
  _Tp
  __comp_ellint_1_agm(_Tp __k)
  {
    constexpr auto _S_max_iter = 100;
    auto __am = (_Tp{1} + __k) / _Tp{2};
    auto __gm = std::sqrt(__k);
    for (int __i = 0; __i < _S_max_iter; ++__i)
      {
	__gm = std::sqrt(__gm * std::exchange(__am, (__am + __gm) / _Tp{2}));
	if (std::abs(__am - __gm) < std::numeric_limits<_Tp>::epsilon())
	  break;
      }
    return __gm;
  }

/**
 *  Use MacLaurin series to calculate the elliptic nome
 *  given the , k.
 */
template<typename _Tp>
  _Tp
  __ellint_nome_series(_Tp __k)
  {
    auto __m = __k * __k; 
    return __m * ((_Tp{1} / _Tp{16})
	 + __m * ((_Tp{1} / _Tp{32})
	 + __m * ((_Tp{21} / _Tp{1024})
	 + __m * ((_Tp{31} / _Tp{2048})
	 + __m * (_Tp{6257} / _Tp{524288})))));
  }

/**
 *  Use the arithmetic-geometric mean to calculate the elliptic nome
 *  given the , k.
 */
template<typename _Tp>
  _Tp
  __ellint_nome_agm(_Tp __k)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(std::real(__k));
    auto __kp = std::sqrt((_Tp{1} - __k) * (_Tp{1} + __k));
    auto __K = __comp_ellint_1_agm(__k);
    auto __Kp = __comp_ellint_1_agm(__kp);
    return std::exp(-_S_pi * __Kp / __K);
  }

/**
 *  Return the elliptic nome given the , k.
 */
template<typename _Tp>
  _Tp
  __ellint_nome(_Tp __k)
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (__k < std::pow(_Tp{67} * _S_eps, _Tp{1}/_Tp{8}))
      return __ellint_nome_series(__k);
    else
      return __ellint_nome_agm(__k);
  }


template<typename _Tp>
  void
  test_K(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));

    auto width = 6 + std::cout.precision();
    const auto del = _Tp{1} / _Tp{1000};
    for (int i = 0; i <= 1000; ++i)
      {
	auto k = i * del;
	std::cout << ' ' << std::setw(width) << k
		  << ' ' << std::setw(width) << __comp_ellint_1_agm(k)
		  << ' ' << std::setw(width) << __ellint_nome(k) << '\n';
      }
  }


int
main()
{
  std::cout << "\n\nfloat\n";
  test_K<float>();

  std::cout << "\n\ndouble\n";
  test_K<double>();

  std::cout << "\n\nlong double\n";
  test_K<long double>();

  std::cout << "\n\n__float128\n";
  test_K<__float128>();
}
