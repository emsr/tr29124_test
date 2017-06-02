// This is a demo program for the jacobi library

/*
$HOME/bin/bin/g++ -o test_gauss_jacobi test_gauss_jacobi.cpp gauss_jacobi_interface.cpp

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_gauss_jacobi test_gauss_jacobi.cpp gauss_jacobi_interface.cpp
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "jacobi.h"

double 
fun(double x)
{
  return std::cos(3.0 * x);
}

int 
main(int argc, char **argv)
{
  int Q;
  
  if (argc == 1)
    {
      printf("Usage: integrate n\nn is the number of quadrature points\n");
      return 1;
    }
  else
    Q = atoi(argv[1]);
  

  if (Q < 1)
    {
      printf("Usage: integrate n\nn is the number of quadrature points\n");
      return 1;
    }

  double *f = (double *)malloc(4 * Q * sizeof(double));
  double *ws = f + Q;

  jac_quadrature<double> *quad = jac_quadrature_alloc<double>(Q);
  jac_quadrature_zwd(quad, Gauss, 0.0, 0.0, ws);

  for (int i = 0; i < Q; ++i)
    f[i] = fun(quad->x[i]);

  double integr = jac_integrate(quad, f);
  double exact = 2.0 / 3.0 * std::sin(3.0);

  printf("Integral of cos(3x) from -1 to 1: %lf\n", integr);
  printf("Error: %e\n", integr - exact);
  jac_quadrature_free(quad);
  free(f);

  return 0;
}
