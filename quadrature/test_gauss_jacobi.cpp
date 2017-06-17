// This is a demo program for the jacobi library

/*
$HOME/bin/bin/g++ -o test_gauss_jacobi test_gauss_jacobi.cpp

$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_gauss_jacobi test_gauss_jacobi.cpp
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "jacobi.h"

double
fun(double x)
{
  return std::cos(3.0 * x);
}

template<typename _Tp>
  void
  test_gauss_jacobi(int Q, _Tp alpha, _Tp beta)
  {
    jac_quadrature<double> quad(Gauss, Q, alpha, beta);

    auto integr = quad.integrate(fun);
    auto exact = 2.0 / 3.0 * std::sin(3.0);

    std::cout << "Integral of cos(3x) from -1 to 1: " << integr << '\n';
    std::cout << "Error: " << integr - exact << '\n';
  }


int 
main(int n_app_args, char **app_arg)
{
  int Q = 10;
  double alpha = 0.0, beta = 0.0;
  if (n_app_args == 1)
    {
      std::cerr << "Usage: integrate n [alpha] [beta]\n"
		<< "n is the number of quadrature points\n"
		<< "alpha and beta are optional Jacobi parameters.\n";
      return 1;
    }
  else
    {
      Q = atoi(app_arg[1]);
      if (Q < 1)
	{
	  printf("Usage: integrate n\nn is the number of quadrature points\n");
	  return 1;
	}
      if (n_app_args == 3)
	beta = alpha = strtod(app_arg[2], nullptr);
      else
	{
	  alpha = strtod(app_arg[2], nullptr);
	  beta = strtod(app_arg[3], nullptr);
	}
    }

  test_gauss_jacobi(Q, alpha, beta);

  return 0;
}
