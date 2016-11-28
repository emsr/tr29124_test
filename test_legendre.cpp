/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include "legendre.tcc"
#include "bits/sf_legendre.tcc"

  template<typename _Tp>
    std::vector<_Tp>
    __legendre_zeros(unsigned int __l, _Tp proto = _Tp{})
    {
      const auto _S_eps = __gnu_cxx::__epsilon(proto);
      const auto _S_pi = __gnu_cxx::__const_pi(proto);

      auto __m = __l / 2;
      std::vector<_Tp> __zero(__l);
      if (__l & 1)
	__zero[__m] = _Tp{0};
      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  // Clever approximation of root.
	  auto __z = std::cos(_S_pi * (__i - _Tp{1} / _Tp{4})
				    / (__l + _Tp{1} / _Tp{2}));
	  auto __k = 0;
	  auto __z1 = __z;
	  do
	    {
	      // Compute __P, __P1, and __P2 the Legendre polynomials of order
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute __Pp the derivative of the Legendre polynomial of order l.
	      auto __P1 = _Tp{0};
	      auto __P = _Tp{1};
	      for  (auto __j = 1u; __j <= __l; ++__j)
		{
		  auto __P2 = __P1;
		  __P1 = __P;
		  // Recursion for legendre polynomials.
		  __P = ((_Tp{2} * __j - _Tp{1}) * __z * __P1
		      - (__j - _Tp{1}) * __P2) / __j;
		}
	      // Recursion for the derivative of The Legendre polynomials.
	      auto __Pp = __l * (__z * __P - __P1) / (__z * __z - _Tp{1});
	      __z1 = __z;
	      // Converge on root by Newton's method.
	      __z = __z1 - __P / __Pp;
	      ++__k;
	    }
	  while (std::abs(__z - __z1) > _S_eps);

	  __zero[__i - 1] = -__z;
	  __zero[__l - __i] = __z;
	  //__w[__i - 1] = _Tp{2} / ((_Tp{1} - __z * __z) * __Pp * __Pp);
	  //__w[__l - __i] = __w[__i - 1];
	}
      return __zero;
    }

template<typename _Tp>
  void
  test_legendre(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int l = 0; l <= 50; ++l)
      {
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << "P_" << l << "(x)";
	std::cout << '\n';
	for (int i = -100; i <= 100; ++i)
	  {
	    double x = i * _Tp{0.01Q};
	    double P, P_l, P_l0;
	    try
	      {
		P = __legendre_p(l, x);
		P_l = std::__detail::__poly_legendre_p(l, x);
		P_l0 = std::__detail::__assoc_legendre_p(l, 0, x);
	      }
	    catch (std::exception & err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    std::cout << "  " << std::setw(width) << x;
	    std::cout << "  " << std::setw(width) << P;
	    std::cout << "  " << std::setw(width) << P_l;
	    std::cout << "  " << std::setw(width) << P_l0;
	    std::cout << "  " << std::setw(width) << P_l - P;
	    std::cout << '\n';
	  }
      }

    for (unsigned int l = 0; l <= 3; ++l)
      {
	for (unsigned int m = l; m <= l; --m)
	  {
	    std::cout << "  " << std::setw(width) << "x";
	    std::cout << "  " << std::setw(width) << "P_" << l << "_" << m << "(x)";
	    std::cout << '\n';
	    for (int i = -100; i <= 100; ++i)
	      {
		double x = i * _Tp{0.01Q};
		double P_lm;
		try
		  {
		    P_lm = std::__detail::__assoc_legendre_p(l, m, x);

		  }
		catch (std::exception & err)
		  {
		    std::cerr << err.what() << '\n';
		  }
		std::cout << "  " << std::setw(width) << x;
		std::cout << "  " << std::setw(width) << P_lm;
		std::cout << '\n';
	      }
	  }
      }

    for (int l = 0; l <= 50; ++l)
      {
	auto zero = __legendre_zeros(l, proto);
	std::cout << "\nl = " << std::setw(4) << l << ":\n";
	for (auto z : zero)
	  std::cout << ' ' << std::setw(width) << z << '\n';
      }
  }

int
main()
{
  test_legendre(1.0);

  return 0;
}

