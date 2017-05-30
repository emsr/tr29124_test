#include"Bernstein.h"

// Testing for the Bernstein library: 
//   Root solver 1st case: multiple roots
//
// by Evan Yi-Feng Tsai, 2001
void main(void)
{
  int ii, nn = 20;
  double error, delta, eta, epsilon;
  Bernstein PP, roots;
  FILE *out;

  // Set the tolerances:
  epsilon = 1e-5;
  delta = 1e-11;
  eta = 1e-8;

  //  Set up the polynomial:
  //    This polynomial has very special roots: 
  //    If the polynomial is of degree n, then 
  //    there are two single roots at 0 and 1,
  //    and also (n-2) multiple roots at 0.5    
  PP.dgr = nn;
  try {  PP.cf = new double[PP.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In file test_rootsolver.cpp: " << exception.what() << endl;
    exit(1);
  }

  for ( ii = 0; ii <= PP.dgr; ii++ )
    PP.cf[ii] = pow(-1.0, ii) * ii * (nn - ii); 

  out = fopen("res8","w");

  roots = ROOT(PP, delta, eta, epsilon);
  if ( roots.dgr == -1 ) {
    cout << "Error occurred in root solving!  Method aborts." << endl;
    exit(1);
  }

  fprintf(out,"Test case -- Rootsolver (first example):\n\n");
  // Verify the roots found:
  for ( ii = 0; ii < PP.dgr; ii++ )
  {
    error = EVAL(PP, roots.cf[ii]);
    fprintf(out, "The error for root %1.16f = %1.16f\n", roots.cf[ii], error);  
  }

  fclose(out);
}
