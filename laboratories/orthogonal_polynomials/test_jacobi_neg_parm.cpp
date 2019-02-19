/*
$HOME/bin/bin/g++ -std=c++2a -g -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -Wall -Wextra -Wno-psabi -o test_jacobi_neg_parm test_jacobi_neg_parm.cpp
./test_jacobi_neg_parm > test_jacobi_neg_parm.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

template<typename Tp>
  void
  test_jacobi_neg_parm(int n, Tp alpha, Tp beta)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    std::cout << '\n' << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * Tp{0.01L};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << __gnu_cxx::jacobi(n, alpha - 0.01, beta, x)
		  << ' ' << std::setw(w) << __gnu_cxx::jacobi(n, alpha, beta, x)
		  << ' ' << std::setw(w) << __gnu_cxx::jacobi(n, alpha + 0.01, beta, x)
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_jacobi_neg_parm_alpha(int n, Tp alpha)
  {
    for (int m = 0; m <= 2; ++m)
      {
	Tp beta = -n - alpha - m;
	test_jacobi_neg_parm(n, alpha, beta);
      }
  }

int
main()
{
  int n = 10;
  double alpha = double(2);
  test_jacobi_neg_parm_alpha(n, alpha);
}
