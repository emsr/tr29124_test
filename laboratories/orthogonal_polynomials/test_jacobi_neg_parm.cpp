/*
$HOME/bin/bin/g++ -std=c++2a -g -I. -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -Wall -Wextra -Wno-psabi -o test_jacobi_neg_parm test_jacobi_neg_parm.cpp
./test_jacobi_neg_parm > test_jacobi_neg_parm.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>
#include <sf_jacobi_neg_params.tcc>

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
		  << ' ' << std::setw(w) << lab::__jacobi_recur(n, alpha - 0.001, beta, x).__P_n
		  << ' ' << std::setw(w) << lab::__jacobi_recur(n, alpha, beta, x).__P_n
		  << ' ' << std::setw(w) << lab::__jacobi_recur(n, alpha + 0.001, beta, x).__P_n
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