#include"Bernstein.h"
#include <time.h>

// Demonstration of constructing a degree n
// Chebyshev polynomial and finding its roots:
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii, nn;
  nn = 30;    // the polynomial degree
  double zero, *zeros, pi = 3.1415926535897932;
  double RR, delta, eta;
  Bernstein *TT, Cheb, result;
  FILE *out;

  clock_t start, finish;
  out = fopen("res1","w");

  // Set the tolerances:
  delta = 1e-11;
  eta = 1e-8;

  try {  
    TT = new Bernstein[nn+1];
    TT[0].cf = new double[1];
    TT[1].cf = new double[2]; 
    zeros = new double[nn];
  }
  catch ( bad_alloc exception ) {  
    cout << "In test_Chebyshev.cpp: " << exception.what() << endl;  
    exit(1);
  }

  start = clock();

  // Set T[0] = 1: 
  TT[0].dgr = 0;
  TT[0].cf[0] = 1.0;

  // Set T[1] = 2u-1 = (-1)(1-u) + 1(u):
  TT[1].dgr = 1;
  TT[1].cf[0] = -1.0; 
  TT[1].cf[1] = 1.0;

  // Recursion for constructing Tn:
  for ( ii = 2; ii <= nn; ii++ )
    TT[ii] = 2.0 * TT[1] * TT[ii-1] - TT[ii-2];

  Cheb = TT[nn]; 

  // Root solving (we know there are only single roots):
  result = sort ( sROOT(Cheb, delta, eta) );
  if ( result.dgr == -1 ) {
    cout << "In test_Chebyshev.cpp: " << "Error in root solving." << endl;
    exit(1);
  }

  // Compute ``standard'' roots:
  for ( ii = nn - 1; ii >= 0; ii-- )
    zeros[nn-1-ii] = 0.5 + 0.5 * 
                   cos((2.0*(double)ii+1.0)*pi/(2.0*(double)nn));

  // Compute RMS error:
  RR = 0.0;
  for ( ii = 0; ii <= nn - 1; ii++ ) {
    zero = 0.5 + 0.5 * 
           cos((2.0*(double)(nn-1-ii)+1.0)*pi/(2.0*(double)nn));
    RR = RR + pow((result.cf[ii] - zero),2.0);
  }
  RR = RR / (double)nn;
  RR = sqrt(RR);  

  finish = clock();

  // Print the results:
  fprintf(out,"Test case -- Root solving for Chebyshev polynomials:\n\n");
  fprintf(out,"Roots of a degree %d Chebyshev polynomial:\n\n", nn);
  fprintf(out,"Computed            Standard            Error\n");
  for ( ii = 0; ii <= nn - 1; ii++ )
    fprintf(out,"%1.16f  %1.16f  %1.16f\n",result.cf[ii],zeros[ii],result.cf[ii]-zeros[ii]);
  
  fprintf(out,"\nRMS error: %1.16f\n", RR);

  fprintf(out,"Time taken: %f seconds\n",((double) (finish-start)) / CLOCKS_PER_SEC);

  fclose(out);
}






