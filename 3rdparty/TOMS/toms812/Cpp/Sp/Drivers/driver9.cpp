#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Root solver 2nd case: arbitrary roots
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii;
  double rr[10];
  double aa[10][2];
  double error, delta, eta, epsilon;
  Bernstein AA[10], xx, roots;
  FILE *out;

  // The roots:
  rr[0] = 0.15; 
  rr[1] = 0.24;
  rr[2] = 0.39;
  rr[3] = 0.48; 
  rr[4] = 0.54;
  rr[5] = 0.66;
  rr[6] = 0.72;
  rr[7] = 0.81; 
  rr[8] = 0.93;
  rr[9] = 0.05; 

  // Set the tolerances:
  epsilon = 1e-5;
  delta = 1e-11;
  eta = 1e-8;

  // Set up the original polynomial:
  for ( ii = 0; ii < 10; ii++ ) {
    aa[ii][0] = -1.0 * rr[ii];
    aa[ii][1] = 1.0 - rr[ii];
    AA[ii] = Bernstein(*(aa+ii),1);
  } 

  xx = AA[1]*AA[2]*AA[3]*AA[4]*AA[5]*AA[6]*AA[7]*AA[8]*AA[9];  
  //  xx = (AA[0]^3)*AA[1]*AA[2]*AA[3]*AA[4]*(AA[5]^2)*AA[6]*AA[7]*(AA[8]^4)*AA[9];

  out = fopen("res9","w");

  roots = ROOT (xx, delta, eta, epsilon);
  if ( roots.dgr == -1 ) {
    cout << "Error occurred in root solving!  Method aborts." << endl;
    exit(1);
  }

  //  Verify the roots found:
  //    The error is determined by plugging in the computed roots
  //    to the polynomial to see how far does the result deviate
  //    from zero.
  fprintf(out,"Test case -- Rootsolver (second example):\n\n");
  for ( ii = 0; ii < roots.dgr; ii++ )
  {
    error = EVAL(xx, roots.cf[ii]);
    fprintf(out,"The error for root %1.16f = %1.16f\n",roots.cf[ii],error);  
  }

  fclose(out);
}
