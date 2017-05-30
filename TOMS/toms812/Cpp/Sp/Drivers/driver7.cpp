#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Normalization
//
// The RMS value of a normalized polynomial on [0,1] should
// have a numerical value of 1 
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int nn = 15;
  double aa[16]; 
  double RR;

  // Set up the original polynomial:
  aa[0] = 6.3; aa[1] = -3.42; aa[2] = 0.06; 
  aa[3] = -3.9; aa[4] = 1.42; aa[5] = 6.5;
  aa[6] = 4.2; aa[7] = -0.42; aa[8] = 0.06;
  aa[9] = 31.3; aa[10] = -1.42; aa[11] = 0.6; 
  aa[12] = 3.9; aa[13] = -10.42; aa[14] = -6.5;
  aa[15] = 40.2;

  Bernstein xx(aa, nn), result, value;

  FILE *out;
  out = fopen("res7","w");

  result = Normalize(xx);

  // Obtain the RMS polynomial:
  value = integrate(result * result);
  RR = value.cf[value.dgr];
  fprintf(out,"Test case -- Normalization:\n\n");
  fprintf(out,"The RMS value of the normalized polynomial on [0,1]: %1.16f\n", RR);
  fprintf(out,"Its error from the theoretical value is: %1.16f\n", RR - 1.0);  
  fclose(out);
}

