#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Subdivision by de Castaljau algorithm
//
// This function verifies the subdivision algorithm
// by comparing points on the subdivided curves with
// the corresponding points on the original curve.
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int nn = 5;
  double xi, xiL, xiR;
  double xil, xir;
  Bernstein rr, left, right;
  FILE *out;

  out = fopen("res10","w");

  rr.dgr = nn;
  try {  rr.cf = new double[rr.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In test_subdivision.cpp: " << exception.what() << endl;
    exit(1);
  }

  rr.cf[0] = 65.8;  rr.cf[1] = -34.56;  rr.cf[2] = 11.004;
  rr.cf[3] = -51.72;  rr.cf[4] = 30.6;  rr.cf[5] = -87.0;

  // Enter the left testing point, the subdivision point,
  // and the right testing point:
  xiL = 0.12; xi = 0.38; xiR = 0.95;

  fprintf(out,"Test case -- Subdivision:\n\n");
  fprintf(out,"The subdivision point: %f\n", xi);
  fprintf(out,"The left test point (in the original domain [0,1]): %f\n", xiL);
  fprintf(out,"The right test point (in the original domain [0,1]): %f\n", xiR);

  left = subLEFT(rr,xi);
  right = subRIGHT(rr,xi);

  // xi values on the local coordinates:
  xil = xiL / xi;
  xir = (xiR - xi) / (1.0 - xi);

  fprintf(out,"\nError on the left piece: %1.16f\n",EVAL(left,xil)-EVAL(rr,xiL));
  fprintf(out,"Error on the right piece: %1.16f\n",EVAL(right,xir)-EVAL(rr,xiR));

  fclose(out);
}
