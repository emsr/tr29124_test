/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_gegenbauer test_gegenbauer.cpp -lquadmath
./test_gegenbauer > test_gegenbauer.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_gegenbauer test_gegenbauer.cpp -lquadmath
./test_gegenbauer > test_gegenbauer.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "jacobi_small.hpp"
#include <bits/float128_io.h>

  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __gegenbauer_zeros(unsigned int __n, _Tp __alpha)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha);
      const unsigned int _S_maxit = 1000u;
      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);

      _Tp __z;
      _Tp __w;
      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  if (__i == 1)
	    {
	      auto __an = __alpha / __n;
	      auto __an2 = __an * __an;
	      auto __r1 = (1.0 + __alpha) * (2.78 / (4.0 + __n * __n)
			+ 0.768 * __an / __n);
	      auto __r2 = 1.0 + 1.48 * __an + 0.96 * __an + 1.282 * __an2;
	      __z = 1.0 - __r1 / __r2;
	    }
	  else if (__i == 2)
	    {
	      auto __r1 = (4.1 + __alpha)
			/ ((1.0 + __alpha) * (1.0 + 0.156 * __alpha));
	      auto __r2 = 1.0
			+ 0.06 * (__n - 8.0) * (1.0 + 0.12 * __alpha) / __n;
	      auto __r3 = 1.0
		    + 0.012 * __alpha * (1.0 + 0.25 * std::abs(__alpha)) / __n;
	      __z -= (1.0 - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == 3)
	    {
	      auto __r1 = (1.67 + 0.28 * __alpha) / (1.0 + 0.37 * __alpha);
	      auto __r2 = 1.0 + 0.22 * (__n - 8.0) / __n;
	      auto __r3 = 1.0 + 8.0 * __alpha / ((6.28 + __alpha) * __n * __n);
	      __z -= (__pt[0].__zero - __z) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n - 1)
	    {
	      auto __r1 = (1.0 + 0.235 * __alpha) / (0.766 + 0.119 * __alpha);
	      auto __r2 = 1.0 / (1.0 + 0.639 * (__n - 4.0)
						/ (1.0 + 0.71 * (__n - 4.0)));
	      auto __r3 = 1.0 / (1.0 + 20.0 * __alpha
				/ ((7.5 + __alpha) * __n * __n));
	      __z += (__z - __pt[__n - 4].__zero) * __r1 * __r2 * __r3;
	    }
	  else if (__i == __n)
	    {
	      auto __r1 = (1.0 + 0.37 * __alpha) / (1.67 + 0.28 * __alpha);
	      auto __r2 = 1.0 / (1.0 + 0.22 * (__n - 8.0) / __n);
	      auto __r3 = 1.0 / (1.0 + 8.0 * __alpha
				 / ((6.28 + __alpha) * __n * __n));
	      __z += (__z - __pt[__n - 3].__zero) * __r1 * __r2 * __r3;
	    }
	  else
	    __z = 3.0 * __pt[__i - 2].__zero
		- 3.0 * __pt[__i - 3].__zero + __pt[__i - 4].__zero;

	  auto __2alpha = _Tp{2} * __alpha;
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __temp = _Tp{2} + __2alpha;
	      auto __C1 = (__temp * __z) / _Tp{2};
	      auto __C2 = _Tp{1};
	      for (auto __j = 2u; __j <= __n; ++__j)
		{
		  auto __C3 = __C2;
		  __C2 = __C1;
		  __temp = _Tp{2} * __j + __2alpha;
		  auto __a = _Tp{2} * __j * (__j + __2alpha)
			   * (__temp - _Tp{2});
		  auto __b = (__temp - _Tp{1})
			   * __temp * (__temp - _Tp{2}) * __z;
		  auto __c = _Tp{2} * (__j - 1 + __alpha)
			   * (__j - 1 + __alpha) * __temp;
		  __C1 = (__b * __C2 - __c * __C3) / __a;
		}
	      auto __Cp = (__n * (-__temp * __z) * __C1
			   + _Tp{2} * (__n + __alpha) * (__n + __alpha) * __C2)
			/ (__temp * (_Tp{1} - __z * __z));
	      auto __z1 = __z;
	      __z = __z1 - __C1 / __Cp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = std::exp(std::lgamma(__alpha + _Tp(__n))
			       + std::lgamma(__alpha + _Tp(__n))
			       - std::lgamma(_Tp(__n + 1))
			       - std::lgamma(_Tp(__n + 1) + __2alpha))
		      * __temp * std::pow(_Tp{2}, __2alpha) / (__Cp * __C2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__jacobi_zeros: Too many iterations");
	    }
	  __pt[__i - 1].__zero = __z;
	  __pt[__i - 1].__weight = __w;
	}

      return __pt;
    }

template<typename _Tp>
  void
  test_gegenbauer(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 6;

    std::cout << "\njacobi\n";
    for (int n = 0; n <= 5; ++n)
      {
	for (int i = 0; i <= 3; ++i)
	  {
            auto alpha = i * _Tp{1.0Q};
            auto jalpha = alpha - _Tp{0.5Q};
            auto jnorm = __gnu_cxx::pochhammer(_Tp{2} * alpha, n)
			/ __gnu_cxx::pochhammer(alpha + _Tp{0.5Q}, n);
            std::cout << "n     = " << n << '\n';
            std::cout << "alpha = " << alpha << '\n';
            Life::Jacobi<_Tp> jac(n, jalpha, jalpha);
	    for (int k = 0; k <= 200; ++k)
              {
        	auto x = (k - 100) * _Tp{0.01Q};
        	auto ggb = __gnu_cxx::gegenbauer(n, alpha, x);
        	auto gj = jnorm * __gnu_cxx::jacobi(n, jalpha, jalpha, x);
        	auto lj = jnorm * jac(x);
        	std::cout << std::setw(width) << x
        	          << std::setw(width) << ggb
        	          << std::setw(width) << gj
        	          << std::setw(width) << lj
        	          << std::setw(width) << ggb - lj
        	          << '\n';
              }
            std::cout << '\n';

	    for (int n = 0; n <= 50; ++n)
	      {
		auto pt = __gegenbauer_zeros(n, alpha);
		std::cout << "\nn = " << std::setw(4) << n
			  << ", alpha = " << std::setw(width) << alpha << ":\n";
		for (auto [z, w] : pt)
		  std::cout << ' ' << std::setw(width) << z
		  	    << ' ' << std::setw(width) << w
		  	    << '\n';
	      }
          }
      }
  }

int
main()
{
  test_gegenbauer(1.0F);

  test_gegenbauer(1.0);

  test_gegenbauer(1.0L);

  test_gegenbauer(1.0Q);

  return 0;
}
