#include"Bernstein.h"
#include <time.h>

// Testing for the Bernstein library: 
//   Composition operation: <<
//
// Caution: The testing of this operation 
// is subject to extrapolation, if the 
// range of y(x) is beyond [0,1]. 
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii, mm = 15, nn = 9;
  double aa[16], bb[10]; 
  double xi, dxi = 1.0 / 100.0, error, RR;
  clock_t start, finish;

  // Set up the two polynomials:
  aa[0] = 6.3; aa[1] = -3.42; aa[2] = 0.06; 
  aa[3] = -3.9; aa[4] = 1.42; aa[5] = 6.5;
  aa[6] = 4.2; aa[7] = -0.42; aa[8] = 0.06;
  aa[9] = 31.3; aa[10] = -1.42; aa[11] = 0.6; 
  aa[12] = 3.9; aa[13] = -10.42; aa[14] = -6.5;
  aa[15] = 40.2;

  bb[0] = 0.5; bb[1] = 0.246; bb[2] = 0.51; 
  bb[3] = -0.69; bb[4] = 0.44; bb[5] = 0.99;
  bb[6] = 0.51; bb[7] = -0.246; bb[8] = 0.32; 
  bb[9] = -0.74;

  Bernstein x(aa, mm), y(bb, nn), result;

  FILE *out;
  out = fopen("res4","w");

  fprintf(out,"Test case -- Composition operation:\n\n");
  fprintf(out,"For the two Bernstein-form polynomials x and y,\n");
  fprintf(out,"where deg(x) = %d, deg(y) = %d\n\n", mm, nn);

  start = clock();
  result = x << y;
  finish = clock();

  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = (EVAL(result,xi)) - EVAL(x,EVAL(y,xi));
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);

  fprintf(out,"The RMS error in computing x << y: %1.16f\n",RR);
  fprintf(out,"Time taken: %f seconds\n",((double) (finish-start)) / CLOCKS_PER_SEC);
  fclose(out);
}
