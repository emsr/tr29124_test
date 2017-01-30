/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I.. -o build_clenshaw_curtis build_clenshaw_curtis.cpp -lquadmath
./build_clenshaw_curtis > build_clenshaw_curtis.txt
*/

#include <vector>
#include <ext/cmath>
#include <bits/specfun_state.h>
#include <iostream>
#include <iomanip>

  /**
   * Return a vector of zeros of the Chebyshev function of the second kind
   * of order @f$ n @f$, @f$ U_n(x) @f$.
   * The zeros are given by:
   * @f[
   *   x_k = \cos\left(\frac{k\pi}{n+1}\right), k \elem {1, ..., n}
   * @f]
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __chebyshev_u_zeros(unsigned int __n)
    {
      const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);
      for (unsigned int __k = 1; __k <= __n; ++__k)
	{
	  auto __arg = _Tp(__k) / _Tp(__n + 1);
	  auto __half = __gnu_cxx::__fp_is_equal<_Tp>(__arg, _Tp{0.5L});
	  auto __z = (__half ? _Tp{0} : __gnu_cxx::cos_pi(__arg));
	  auto __w = _S_pi * (_Tp{1} - __z * __z) / _Tp(__n + 1);
	  __pt[__k - 1].__zero = __z;
	  __pt[__k - 1].__weight = __w;
	}
      return __pt;
    }

template<typename _Tp, size_t _N>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_clenshaw_curtis()
  {
    const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
    auto uz = __chebyshev_u_zeros<_Tp>(_N - 1);
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_N + 1);
    __out[0].__zero = _Tp{-1};
    {
      auto __sum = _Tp{0};
      for (auto __j = 1u; __j <= _N / 2; ++__j)
	{
	  auto __b = _Tp(__j == _N / 2 ? 1 : 2);
	  __sum += __b / _Tp(4 * __j * __j - 1);
	}
      __out[0].__weight = (_Tp{1} - __sum) / _Tp{_N};
    }
    for (auto __k = 1u; __k <= uz.size(); ++__k)
      {
	__out[__k].__zero = uz[__k - 1].__zero;

	auto __sum = _Tp{0};
	for (auto __j = 1u; __j <= _N / 2; ++__j)
	  {
	    auto __b = _Tp(__j == _N / 2 ? 1 : 2);
	    __sum += __b * std::cos(2 * __j * __k * _S_pi / _N)
		   / _Tp(4 * __j * __j - 1);
	  }
	auto __w = _Tp{2} * (_Tp{1} - __sum) / _Tp{_N};
	__out[__k].__weight = __w;
      }
    __out[_N].__zero = _Tp{1};
    __out[_N].__weight = __out[0].__weight;
    return __out;
  }

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10<long double>());
  auto w = 8 + std::cout.precision();

  auto cc24 = build_clenshaw_curtis<long double, 24>();
  std::cout << '\n';
  for (const auto& cc : cc24)
    std::cout << std::setw(w) << cc.__zero << ' '
	      << std::setw(w) << cc.__weight << '\n';

  std::cout << '\n';
  auto cc48 = build_clenshaw_curtis<long double, 48>();
  for (const auto& cc : cc48)
    std::cout << std::setw(w) << cc.__zero << ' '
	      << std::setw(w) << cc.__weight << '\n';
}
