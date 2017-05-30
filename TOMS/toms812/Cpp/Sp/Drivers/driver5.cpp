#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Differentiation & integration
//
// This function first integrates a polynomial x,
// and then differentiate it.  The result should
// be identical to x.
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii, nn = 15;
  double aa[16]; 
  double xi, dxi = 1.0 / 100.0, error, RR;

  // Set up the original polynomial:
  aa[0] = 6.3; aa[1] = -3.42; aa[2] = 0.06; 
  aa[3] = -3.9; aa[4] = 1.42; aa[5] = 6.5;
  aa[6] = 4.2; aa[7] = -0.42; aa[8] = 0.06;
  aa[9] = 31.3; aa[10] = -1.42; aa[11] = 0.6; 
  aa[12] = 3.9; aa[13] = -10.42; aa[14] = -6.5;
  aa[15] = 40.2;

  Bernstein xx(aa,nn), int_diff;

  FILE *out;
  out = fopen("res5","w");

  int_diff = diff(integrate(xx));

  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(int_diff,xi) - EVAL(xx,xi); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);

  fprintf(out,"Test case -- Differentiation and integration:\n\n");
  fprintf(out,"The RMS error in integrate_differentiate: %1.16f\n",RR);  

  fclose(out);
}
