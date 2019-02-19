/*
$HOME/bin/bin/g++ -std=c++2a -g -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -Wall -Wextra -Wno-psabi -o test_gegenbauer_neg_parm test_gegenbauer_neg_parm.cpp
./test_gegenbauer_neg_parm > test_gegenbauer_neg_parm.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>


template<typename Tp>
  void
  test_gegenbauer_neg_parm(int n, Tp lambda)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    std::cout << '\n' << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * Tp{0.01L};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << __gnu_cxx::gegenbauer(n, lambda - 0.01, x)
		  << ' ' << std::setw(w) << __gnu_cxx::gegenbauer(n, lambda, x)
		  << ' ' << std::setw(w) << __gnu_cxx::gegenbauer(n, lambda + 0.01, x)
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_gegenbauer_neg_parm_lambda(int n)
  {
    for (int m = 0; m <= 2; ++m)
      {
	const auto lambda = Tp(1 - n - m) / Tp{2};
	test_gegenbauer_neg_parm(n, lambda);
      }
  }

int
main()
{
  int n = 10;
  test_gegenbauer_neg_parm_lambda<double>(n);
}
