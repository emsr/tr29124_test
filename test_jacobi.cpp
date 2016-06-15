/*
$HOME/bin_tr29124/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I/usr/local/include/boost -o test_jacobi test_jacobi.cpp
./test_jacobi > test_jacobi.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "jacobi_small.hpp"

template<typename Tp>
  void
  test_jacobi()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto width = std::cout.precision() + 6;
    std::cout << "\njacobi\n";
    for (int n = 0; n <= 5; ++n)
      {
	for (int i = 0; i <= 3; ++i)
	  {
            auto alpha = i * Tp{1.0L};
            for (int j = 0; j <= 3; ++j)
              {
        	auto beta = j * Tp{1.0L};
        	std::cout << "n     = " << n << '\n';
        	std::cout << "alpha = " << alpha << '\n';
        	std::cout << "beta  = " << beta << '\n';
                Life::Jacobi<Tp> jac(n, alpha, beta);
		for (int k = 0; k <= 200; ++k)
        	  {
        	    auto x = (k - 100) * Tp{0.01L};
        	    std::cout << std::setw(width) << x
        	              << std::setw(width) << __gnu_cxx::jacobi(n, alpha, beta, x)
        	              << std::setw(width) << jac(x)
        	              << '\n';
        	  }
        	std::cout << '\n';
	      }
            std::cout << '\n';
          }
      }
  }

int
main()
{
  test_jacobi<long double>();
}
