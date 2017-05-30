#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Arithmetic operations: +, -, *, /, ^
//
// by Evan Yi-Feng Tsai, 2001
// 
void main(void)
{
  int ii, mm = 15, nn = 9;
  double aa[16], bb[10]; 
  double xi, dxi = 1.0 / 100.0, error, RR;

  // Set up the two polynomials:
  aa[0] = 6.3; aa[1] = -3.42; aa[2] = 0.06; 
  aa[3] = -3.9; aa[4] = 1.42; aa[5] = 6.5;
  aa[6] = 4.2; aa[7] = -0.42; aa[8] = 0.06;
  aa[9] = 31.3; aa[10] = -1.42; aa[11] = 0.6; 
  aa[12] = 3.9; aa[13] = -10.42; aa[14] = -6.5;
  aa[15] = 40.2;

  bb[0] = 10.5; bb[1] = 2.46; bb[2] = -5.11; 
  bb[3] = -6.90; bb[4] = 10.0; bb[5] = 1.0;
  bb[6] = 12.5; bb[7] = -2.46; bb[8] = 5.11; 
  bb[9] = -1.90;

  Bernstein xx(aa, mm), yy(bb, nn), result, qq, rr, recovered_xx;

  FILE *out;
  out = fopen("res2","w");

  fprintf(out,"Test case -- Arithmetic operations:\n\n");
  fprintf(out,"For the two Bernstein-form polynomials x and y,\n");
  fprintf(out,"where deg(x) = %d, deg(y) = %d\n\n", mm, nn);

  // Test for the '+' operator:
  result = xx + yy;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - ( EVAL(xx, xi) + EVAL(yy, xi) ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing x + y: %1.16f\n", RR);

  // Test for the '-' operator:
  result = xx - yy;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - ( EVAL(xx, xi) - EVAL(yy, xi) ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing x - y: %1.16f\n", RR);

  // Test for the '*' operator:
  result = xx * yy;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - ( EVAL(xx, xi) * EVAL(yy, xi) ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing x * y: %1.16f\n",RR);

  // Test for the '*' operator --- pre-multiply:
  result = 3.0 * yy;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - ( 3.0 * EVAL(yy, xi) ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing 3.0 * y: %1.16f\n", RR);

  // Test for the '*' operator --- post-multiply:
  result = xx * 3.0;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - ( EVAL(xx, xi) * 3.0 ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing x * 3.0: %1.16f\n", RR);

  // Test for the '/' operator:
  result = xx / yy;
  qq = quo(result, mm-nn);
  rr = rem(result, nn-1);

  // Putting it back:
  recovered_xx = qq * yy + rr;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(recovered_xx, xi) - EVAL(xx, xi); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing x / y: %1.16f\n", RR);  

  // Test for the '^' operator:
  result = yy^2;
  RR = 0.0;
  for ( ii = 0; ii <= 100; ii++ ) {
    xi = (double)ii * dxi;
    error = EVAL(result, xi) - pow( EVAL(yy, xi), 2.0 ); 
    RR += (error * error);
  }

  RR = RR / 101.0;
  RR = sqrt(RR);
  fprintf(out,"The RMS error in computing y^2: %1.16f\n", RR);

  fclose(out);
}





