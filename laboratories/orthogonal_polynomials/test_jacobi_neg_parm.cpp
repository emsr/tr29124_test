/**
 *
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
		  << ' ' << std::setw(w) << lab::jacobi_recur(n, alpha - 0.001, beta, x).P_n
		  << ' ' << std::setw(w) << lab::jacobi_recur(n, alpha, beta, x).P_n
		  << ' ' << std::setw(w) << lab::jacobi_recur(n, alpha + 0.001, beta, x).P_n
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
