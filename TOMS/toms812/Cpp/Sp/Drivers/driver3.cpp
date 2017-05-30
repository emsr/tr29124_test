#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Computation of binomial coefficients
//
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii, nn;
  double aa[2];
  double error, RR;
  Bernstein xx, PP;
  FILE *out;

  out = fopen("res3","w");

  // Creates the trivial polynomial x(t) = 1.0(1-t) + 1.0t:
  aa[0] = 1.0;  aa[1] = 1.0;
  xx = Bernstein(aa, 1);

  nn = 100;
  // Take it to the nth power: all of the Bernstein coefficients
  // of P(t) should have the value 1
  PP = xx^nn;

  fprintf(out,"Test case -- Computation of binomial coefficients:\n\n");
  fprintf(out,"For x(t) = 1.0(1-t) + 1.0t and P(t) = x^%d:\n\n",nn);

  RR = 0.0; 
  for ( ii = 0; ii <= nn; ii++ ) {
    error = PP.cf[ii] - 1.0;
    RR += (error * error);
  }

  RR = RR / (double) (nn+1);
  RR = sqrt(RR);

  fprintf(out,"The RMS error in the Bernstein coefficients of P(t)\n");
  fprintf(out,"(from the theoretical value 1) = %1.18f\n",RR);  

  fclose(out);
}
