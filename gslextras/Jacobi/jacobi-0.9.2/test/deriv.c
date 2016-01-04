// Function to calculate the derivative. It will output the maximum error of the derivatives

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/jacobi.h"



double 
fun(double x)
{
  return cos(3.0 * x);
}

double 
dfun(double x)
{
  return -3.0*sin(3.0*x);
}


int 
main(int argc, char **argv)
{
  int Q;
  
  if (argc == 1)
    {
      printf("Usage: deriv n\nn is the number of quadrature points\n");
      return 1;
    }
  else
    Q = atoi(argv[1]);
  

  if (Q < 1)
    {
      printf("Usage: integrate n\nn is the number of quadrature points\n");
      return 1;
    }

  double *f = (double *) malloc(5*Q * sizeof(double));
  double *d = f + Q;
  
  double *ws = f + 2*Q;

  jac_quadrature *quad = jac_quadrature_alloc(Q);
  jac_quadrature_zwd(quad, JAC_GJ, 0.0, 0.0, ws);
  
  
  int i;

  for (i = 0; i < Q; ++i)
    f[i] = fun(quad->x[i]);
  
  jac_differentiate(quad, f, d);

  printf("X \t Derivative \t Error\n");

  for (i = 0; i < quad->Q; ++i)
    printf("%lf \t %lf \t %e\n", quad->x[i], d[i], d[i] - dfun(quad->x[i]));
  
  jac_quadrature_free(quad);
  free(f);
  
  return 0;
}


  
      
