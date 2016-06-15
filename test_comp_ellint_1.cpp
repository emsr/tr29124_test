/*
$HOME/bin_tr29124/bin/g++ -std=gnu++1z -o test_comp_ellint_1 test_comp_ellint_1.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_comp_ellint_1 > test_comp_ellint_1.txt
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <experimental/utility>
#include "float128.h"

//  Use AGM to do an ab initio calculation of K(k).
template<typename _Tp>
  _Tp
  __comp_ellint_1_agm(_Tp __k)
  {
    constexpr auto _S_max_iter = 100;
    auto __am = _Tp{0.5Q} * (_Tp{1} + __k);
    auto __gm = std::sqrt(__k);
    for (int __i = 0; __i < _S_max_iter; ++__i)
      {
	__gm = std::sqrt(__gm * std::exchange(__am, _Tp{0.5Q} * (__am + __gm)));
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
    constexpr auto _S_pi = _Tp{3.1415926535897932384626433832795029Q};
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
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    if (__k < std::pow(_Tp{67} * _S_eps, _Tp{0.125Q}))
      return __ellint_nome_series(__k);
    else
      return __ellint_nome_agm(__k);
  }

template<typename _Tp>
  void
  test_K()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);

    auto width = 6 + std::cout.precision();
    for (int i = 0; i <= 1000; ++i)
      {
	auto k = _Tp(i * 0.001Q);
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
