// Function to interpolate a function. It will output the interpolation and error

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../src/jacobi.h"



double 
fun(double x)
{
  return cos(3.0 * x);
}



int 
main(int argc, char **argv)
{
  int Q;

  double xp[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  int np = 9;
  double fout[9];
  
  
  
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

  double *f = (double *) malloc(4*Q * sizeof(double));
  
  double *ws = f + Q;

  jac_quadrature *quad = jac_quadrature_alloc(Q);
  jac_quadrature_zwd(quad, JAC_GJ, 0.0, 0.0, ws);
  jac_interpmat_alloc(quad, np, xp);

  
  
  int i;

  for (i = 0; i < Q; ++i)
    f[i] = fun(quad->x[i]);
  
  jac_interpolate(quad, f, fout);

  printf("X \t f(x) \t Error\n");

  for (i = 0; i < np; ++i)
    printf("%lf \t %lf \t %e\n", xp[i], fout[i], fout[i] - fun(xp[i]));
  
  
  jac_quadrature_free(quad);
  free(f);
  
  return 0;  

}


  
      
