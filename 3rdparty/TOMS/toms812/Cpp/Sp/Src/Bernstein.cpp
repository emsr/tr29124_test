#include"Bernstein.h"
// Function definitions for BPOLY:
//
// by Evan Yi-Feng Tsai, 2001
//

// ---
// Constructor:
Bernstein::Bernstein()   
{ 
  dgr = 0;
  cf = NULL;
}

// ---
// Constructor with coefficient and degree inputs:
Bernstein::Bernstein(double *coeff, int degree)
{ 
  int kk;
  dgr = degree;

  try {  cf = new double[dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In constructing a Bernstein object: " << exception.what() << endl;  
    dgr = -1;
  }

  for ( kk = 0; kk <= dgr; kk++ )
    cf[kk] = coeff[kk];
}

// ---
// Copy constructor:
Bernstein::Bernstein(const Bernstein &uu)
{
  int kk;
  dgr = uu.dgr;

  try {  cf = new double[dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In constructing a Bernstein object: " << exception.what() << endl;  
    dgr = -1;
  }

  for ( kk = 0; kk <= dgr; kk++ )
    cf[kk] = uu.cf[kk];
}

// ---
// Destructor:
Bernstein::~Bernstein()
{  if ( dgr > -1 )  delete [] cf;  }

// ---
// Function for ``retrieving'' from the ``prime database''.
//   mode 1 : when the client provides information
//     about `n', the degree of the binomial coefficient.
//   mode 2 : when the client provides information
//     about `prNum', the number of primes needed.
//
// This funstion returns a Bernstein instance `output':
//   output.dgr = number of prime numbers,
//   output.cf[] = the prime numbers, in ascending order,
//   output.cf[0] is always 1, not counted as a prime.

Bernstein getprimes(int nn, int prNum, int mode)
{
  int jj;
  Bernstein output;
  double primes[169] =
  {  1.0,  2.0,   3.0,   5.0,   7.0,  11.0,  13.0,  17.0,  19.0,  23.0,  29.0,
    31.0,  37.0,  41.0,  43.0,  47.0,  53.0,  59.0,  61.0,  67.0,  71.0,
    73.0,  79.0,  83.0,  89.0,  97.0, 101.0, 103.0, 107.0, 109.0, 113.0,
   127.0, 131.0, 137.0, 139.0, 149.0, 151.0, 157.0, 163.0, 167.0, 173.0,
   179.0, 181.0, 191.0, 193.0, 197.0, 199.0, 211.0, 223.0, 227.0, 229.0,
   233.0, 239.0, 241.0, 251.0, 257.0, 263.0, 269.0, 271.0, 277.0, 281.0,
   283.0, 293.0, 307.0, 311.0, 313.0, 317.0, 331.0, 337.0, 347.0, 349.0,
   353.0, 359.0, 367.0, 373.0, 379.0, 383.0, 389.0, 397.0, 401.0, 409.0,
   419.0, 421.0, 431.0, 433.0, 439.0, 443.0, 449.0, 457.0, 461.0, 463.0,
   467.0, 479.0, 487.0, 491.0, 499.0, 503.0, 509.0, 521.0, 523.0, 541.0,
   547.0, 557.0, 563.0, 569.0, 571.0, 577.0, 587.0, 593.0, 599.0, 601.0,
   607.0, 613.0, 617.0, 619.0, 631.0, 641.0, 643.0, 647.0, 653.0, 659.0,
   661.0, 673.0, 677.0, 683.0, 691.0, 701.0, 709.0, 719.0, 727.0, 733.0,
   739.0, 743.0, 751.0, 757.0, 761.0, 769.0, 773.0, 787.0, 797.0, 809.0,
   811.0, 821.0, 823.0, 827.0, 829.0, 839.0, 853.0, 857.0, 859.0, 863.0,
   877.0, 881.0, 883.0, 887.0, 907.0, 911.0, 919.0, 929.0, 937.0, 941.0,
   947.0, 953.0, 967.0, 971.0, 977.0, 983.0, 991.0, 997.0  };

  if ( mode == 1 ) {
    prNum = 0;
    for ( jj = 1; jj <= 168; jj++ ) {
      if ( primes[jj] <= nn )  prNum += 1;
      else  break;  
    } 
  }

  output.dgr = prNum;
  try {  output.cf = new double[output.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function getprimes(): " << exception.what() << endl;
    output.dgr = -1;
    return output;
  }

  for ( jj = 0; jj <= output.dgr; jj++ )
    output.cf[jj] = primes[jj];

  return output;
}

// ---
// Function for factorization of a binomial coefficient:
// It returns the exponent array and also the size of this
// array, in an instance of Bernstein object:
Bernstein factorization(int nn, int kk)
{
  int ii, jj;
  double p_i;           
  Bernstein primes, exp;  

  // Irrational case:
  if ( (nn < kk) || (kk < 0) )
    exp.dgr = -1;
  else {
    // Retrieve prime numbers in mode 1
    primes = getprimes(nn, 0, 1);
 
    // Set up the array for storing the exponent vector:
    exp.dgr = primes.dgr;
    try {  exp.cf = new double[exp.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function factorization(): " << exception.what() << endl;
      exp.dgr = -1;
      return exp;
    }

    exp.cf[0] = 1.0;

    for ( jj = 1; jj <= exp.dgr; jj++ ) {
      exp.cf[jj] = 0.0;
      ii = 1; 
      p_i = pow(primes.cf[jj], ii);
      while ( p_i <= nn ) {
        if ( fmod(nn, p_i) < fmod(kk, p_i) )
          exp.cf[jj] += 1.0;
        ii ++; 
        p_i = pow(primes.cf[jj], ii);
      }
    }
  }

  return exp;
}

// ---
// Function for direct evaluation of a binomial coefficient 
// with (n,k) input.
double binom(int nn, int kk)
{
  int jj;
  double result = 1.0; 
  Bernstein primes, exp;

  // Retrieve the prime numbers in mode 1
  primes = getprimes(nn, 0, 1);  

  // Compute the exponents:
  exp = factorization(nn, kk);   

  for ( jj = 1; jj <= exp.dgr; jj++ ) {
    if ( exp.cf[jj] != 0.0 )      
      result *= pow(primes.cf[jj], exp.cf[jj]);  
  }

  return result;
}

// ---
// Function for direct evaluation of a binomial coefficient
// with exponent array input:
double binom(const Bernstein &ee)
{
  int jj;
  double result = 1.0;
  Bernstein primes; 

  // Retrieve the prime numbers in mode 2
  primes = getprimes(0, ee.dgr, 2);

  for ( jj = 1; jj <= ee.dgr; jj++ ) {
    if ( ee.cf[jj] != 0.0 )  
      result *= pow(primes.cf[jj], ee.cf[jj]);
  }

  return result;
}

// ---
// Function for multiplying two binomial coefficients:
// (the output is another exponent array)
Bernstein binomult(const Bernstein &e1, const Bernstein &e2)
{
  int ii, n1, n2;
  Bernstein output;

  n1 = e1.dgr;
  n2 = e2.dgr;

  if ( n1 == -1 || n2 == -1 )
    output.dgr = -1;

  else if ( n1 > n2 ) {
    output.dgr = n1;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binomult(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }

    for ( ii = 1; ii <= n2; ii++ )        
      output.cf[ii] = e1.cf[ii] + e2.cf[ii];  
    for ( ii = n2 + 1; ii <= n1; ii++ )     
      output.cf[ii] = e1.cf[ii];
  }  
  else if ( n1 < n2 ) {
    output.dgr = n2;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binomult(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }

    for ( ii = 1; ii <= n1; ii++ )        
      output.cf[ii] = e2.cf[ii] + e1.cf[ii];  
    for ( ii = n1 + 1; ii <= n2; ii++ )     
      output.cf[ii] = e2.cf[ii];
  }
  else if ( n1 == n2 ) {
    output.dgr = n1;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binomult(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }

    for ( ii = 1; ii <= output.dgr; ii++ )   
      output.cf[ii] = e1.cf[ii] + e2.cf[ii];
  } 
  return output;
}

// ---
// Function for dividing two binomial coefficients:
// (the output is another exponent array)
Bernstein binodiv(const Bernstein &e1, const Bernstein &e2)
{
  int ii, n1, n2;
  Bernstein output;

  n1 = e1.dgr;
  n2 = e2.dgr;

  if ( n1 == -1 || n2 == -1 )
    output.dgr = -1;

  else if ( n1 > n2 )
  {
    output.dgr = n1;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binodiv(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }
    for ( ii = 1; ii <= n2; ii++ )        
      output.cf[ii] = e1.cf[ii] - e2.cf[ii];  
    for ( ii = n2 + 1; ii <= n1; ii++ )     
      output.cf[ii] = e1.cf[ii];
  }  
  else if ( n1 < n2 ) {
    output.dgr = n2;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binodiv(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }
    for ( ii = 1; ii <= n1; ii++ )        
      output.cf[ii] = e1.cf[ii] - e2.cf[ii];  
    for ( ii = n1 + 1; ii <= n2; ii++ )     
      output.cf[ii] = (-1.0) * e2.cf[ii];
  }
  else if ( n1 == n2 ) {
    output.dgr = n1;
    try {  output.cf = new double[output.dgr + 1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function binodiv(): " << exception.what() << endl;
      output.dgr = -1;
      return output;
    }
    for ( ii = 1; ii <= output.dgr; ii++ )   
      output.cf[ii] = e1.cf[ii] - e2.cf[ii];
  }
  return output;
}

// ---
// Function for degree elevation: 
// In the computation for the coefficients, the common
// factor in each summation is factored out to reduce 
// the total number of floating-point multiplications.
Bernstein DE(const Bernstein &uu, int rr)
{
  int jj, jmin, jmax, nn, kk, ii;
  double fracBINOM;
  Bernstein result, e1, e2, e12, e3, e123, CF;

  nn = uu.dgr;
  result.dgr = nn + rr;
  try {  result.cf = new double[result.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function DE(): " << exception.what() << endl;
    result.dgr = -1; 
    return result;
  }
  
  for ( kk = 0; kk <= nn + rr; kk++ ) {
    jmin = kk - rr; 
    if ( jmin < 0 )   jmin = 0;  
    jmax = kk;   
    if ( jmax > nn )   jmax = nn;  
    result.cf[kk] = 0.0;

    // Compute the common factor:
    e1 = factorization(rr, kk-jmin);
    e2 = factorization(nn, jmin);
    e3 = factorization(nn+rr, kk); 
    e12 = binomult(e1, e2);
    e123 = binodiv(e12, e3);
    CF = e123;

    for ( jj = jmin + 1; jj <= jmax; jj++ ) {   
      e1 = factorization(rr, kk-jj);
      e2 = factorization(nn, jj);
      e3 = factorization(nn+rr, kk); 
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);

      for ( ii = 1; ii <= e123.dgr; ii++ ) {
        if ( e123.cf[ii] != CF.cf[ii] )
          CF.cf[ii] = 0.0;
      }
    }

    // Bring out the common factor:
    for ( jj = jmin; jj <= jmax; jj++ ) {    
      e1 = factorization(rr, kk-jj);
      e2 = factorization(nn, jj);
      e3 = factorization(nn+rr, kk); 
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);

      for ( ii = 1; ii <= e123.dgr; ii++ ) {
        if ( e123.cf[ii] == CF.cf[ii] )
          e123.cf[ii] = 0.0;
      }

      fracBINOM = binom(e123);   
      result.cf[kk] += fracBINOM * uu.cf[jj];
    }

    // Multiply the comon factor back:
    result.cf[kk] *= binom(CF);  
  }
  return result;
}

// ---
// Function for evaluating a Bernstein polynomial at x,
// using de Casteljau algorithm:
double EVAL(const Bernstein &uu, double xx)
{
  int nn, kk, rr;
  double result, **tri;

  nn = uu.dgr;

  try { 
    tri = new double *[nn+1];
    for( kk = 0; kk <= nn; kk++ )
      *(tri + kk) = new double[nn+1];  
  }
  catch ( bad_alloc exception ) {  
    cout << "In function EVAL(): " << exception.what() << endl;
    exit(1);
  }

  for ( kk = 0; kk <= nn; kk++ ) 
    tri[kk][0] = uu.cf[kk];  
  for ( rr = 1; rr <= nn; rr++ )
    for ( kk = rr; kk <= nn; kk++ )
      tri[kk][rr] = (1.0 - xx) * tri[kk-1][rr-1] + xx * tri[kk][rr-1];  

  result = tri[nn][nn];

  for ( kk = 0; kk <= nn; kk++ )
    delete [] *(tri + kk);
  delete [] tri;
  return result;
}

// ---
// Function for subdivision of a Bernstein polynomial
// at x and returns the left half:
Bernstein subLEFT(const Bernstein &uu, double xx)
{
  int nn, kk, rr;
  double **tri; 
  Bernstein LEFT;

  nn = uu.dgr; 
 
  try { 
    tri = new double *[nn+1];
    for( kk = 0; kk <= nn; kk++ )
      *(tri + kk) = new double[nn+1];  
  }
  catch ( bad_alloc exception ) {  
    cout << "In function subLEFT(): " << exception.what() << endl;
    LEFT.dgr = -1;
    return LEFT;
  }

  LEFT.dgr = nn;
  try {  LEFT.cf = new double[nn + 1]; }
  catch ( bad_alloc exception ) {  
    cout << "In function subLEFT(): " << exception.what() << endl;
    LEFT.dgr = -1;
    return LEFT;
  }

  for ( kk = 0; kk <= nn; kk++ ) 
    tri[kk][0] = uu.cf[kk];  
  for ( rr = 1; rr <= nn; rr++ )
    for ( kk = rr; kk <= nn; kk++ )
      tri[kk][rr] = (1.0 - xx) * tri[kk-1][rr-1] + xx * tri[kk][rr-1];  

  for ( kk = 0; kk <= nn; kk++ )  
    LEFT.cf[kk] = tri[kk][kk];

  for ( kk = 0; kk <= nn; kk++ )
    delete [] *(tri+kk);
  delete [] tri;
  return LEFT;
}

// ---
// Function for subdivision of a Bernstein polynomial
// at x and returns the right half:
Bernstein subRIGHT(const Bernstein &uu, double xx)
{
  int nn, kk, rr;
  double **tri; 
  Bernstein RIGHT;

  nn = uu.dgr;

  try { 
    tri = new double *[nn+1];
    for( kk = 0; kk <= nn; kk++ )
      *(tri + kk) = new double[nn+1];  
  }
  catch ( bad_alloc exception ) {  
    cout << "In function subLEFT(): " << exception.what() << endl;
    RIGHT.dgr = -1;
    return RIGHT;
  }

  RIGHT.dgr = nn;
  try {  RIGHT.cf = new double[nn + 1]; }
  catch ( bad_alloc exception ) {  
    cout << "In function subLEFT(): " << exception.what() << endl;
    RIGHT.dgr = -1;
    return RIGHT;
  }

  for ( kk = 0; kk <= nn; kk++ ) 
    tri[kk][0] = uu.cf[kk];  
  for ( rr = 1; rr <= nn; rr++ )
    for ( kk = rr; kk <= nn; kk++ )
      tri[kk][rr] = (1.0 - xx) * tri[kk-1][rr-1] + xx * tri[kk][rr-1];  

  for ( kk = 0; kk <= nn; kk++ )  
    RIGHT.cf[kk] = tri[nn][nn-kk];

  for ( kk = 0; kk <= nn; kk++ )
    delete [] *(tri + kk);
  delete [] tri;
  return RIGHT;
}

// ---
// Bernstein "assign" operator:
const Bernstein& Bernstein::operator=(const Bernstein &uu) 
{
  int kk;

  dgr = uu.dgr;

  try {  cf = new double[dgr + 1]; }
  catch ( bad_alloc exception ) {  
    cout << "In function assign(): " << exception.what() << endl;
    dgr = -1;
    return *this;
  }

  for ( kk = 0; kk <= dgr; kk++ )
    cf[kk] = uu.cf[kk];

  return *this;
} 

// ---  
// Bernstein "add" operator:
Bernstein Bernstein::operator+(const Bernstein &uu) const 
{
  int kk, mm, nn;
  mm = dgr;  nn= uu.dgr;
  Bernstein result, temp;
  
  if ( mm == nn ) {
    result.dgr = mm;  
    try {  result.cf = new double[mm+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function add(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    for ( kk = 0; kk <= mm; kk++ )   
      result.cf[kk] = cf[kk] + uu.cf[kk];    
  }
  else if ( mm > nn ) {
    result.dgr = mm;
    try {  result.cf = new double[mm+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function add(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    temp = DE(uu, mm-nn);
    for ( kk = 0; kk <= temp.dgr; kk++ )   
      result.cf[kk] = cf[kk] + temp.cf[kk];  
  }
  else if ( mm < nn ) {
    result.dgr = nn;
    try {  result.cf = new double[nn+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function add(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    temp = DE(*this, nn-mm);
    for ( kk = 0; kk <= temp.dgr; kk++ )   
      result.cf[kk] = temp.cf[kk] + uu.cf[kk];  
  } 
  return result;
} 

// ---  
// Bernstein "subtract" operator:
Bernstein Bernstein::operator-(const Bernstein &uu) const 
{
  int kk, mm, nn;
  mm = dgr;  nn= uu.dgr;
  Bernstein result, temp;
  
  if ( mm == nn ) {
    result.dgr = mm; 
    try {  result.cf = new double[mm+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function subtract(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    for ( kk = 0; kk <= mm; kk++ )   
      result.cf[kk] = cf[kk] - uu.cf[kk];    
  }
  else if ( mm > nn ) {
    result.dgr = mm;
    result.cf = new double[mm+1];
    try {  result.cf = new double[mm+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function subtract(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    temp = DE(uu, mm-nn);
    for ( kk = 0; kk <= temp.dgr; kk++ )   
      result.cf[kk] = cf[kk] - temp.cf[kk];  
  }
  else if ( mm < nn )
  {
    result.dgr = nn;
    try {  result.cf = new double[nn+1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function subtract(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    temp = DE(*this, nn-mm);
    for ( kk = 0; kk <= temp.dgr; kk++ )   
      result.cf[kk] = temp.cf[kk] - uu.cf[kk];  
  } 

  return result;
} 

// ---
// Bernstein "multiply" operator.
//
// In the computation for the coefficients, the common
// factor in each summation is factored out to reduce 
// the total number of floating-point multiplications.
Bernstein Bernstein::operator*(const Bernstein &uu) const 
{
  int jj, jmin, jmax, kk, mm, nn, ii; 
  mm = dgr;  nn= uu.dgr;
  double fracBINOM;
  Bernstein result, e1, e2, e12, e3, e123, CF;

  result.dgr = mm + nn;
  try {  result.cf = new double[result.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function multiply(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }
   
  for ( kk = 0; kk <= mm + nn; kk++ ) {
    jmin = kk - nn;  
    if ( jmin < 0 )    jmin = 0;  
    jmax = kk;      
    if ( jmax > mm )    jmax = mm;  
    result.cf[kk] = 0.0;
  
    // Compute the common factor:
    e1 = factorization(mm, jmin);
    e2 = factorization(nn, kk - jmin);
    e3 = factorization(mm+nn, kk);
    e12 = binomult(e1, e2);
    e123 = binodiv(e12, e3);
    CF = e123;

    for ( jj = jmin + 1; jj <= jmax; jj++ ) {
      e1 = factorization(mm, jj);
      e2 = factorization(nn, kk-jj);
      e3 = factorization(mm+nn, kk); 
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);

      for ( ii = 1; ii <= e123.dgr; ii++ ) {
        if ( e123.cf[ii] != CF.cf[ii] )
        CF.cf[ii] = 0;
      }
    }

    // Bring out the common factor:
    for ( jj = jmin; jj <= jmax; jj++ ) {   
      e1 = factorization(mm, jj);
      e2 = factorization(nn, kk-jj);
      e3 = factorization(mm+nn, kk);
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);

      for ( ii = 1; ii <= e123.dgr; ii++ ) {
        if ( e123.cf[ii] == CF.cf[ii] )
          e123.cf[ii] = 0;
      }

      fracBINOM = binom(e123);   
      result.cf[kk] += fracBINOM * cf[jj] * uu.cf[kk-jj];
    }

    // Multiply the comon factor back:
    result.cf[kk] *= binom(CF);  
  }
  return result;
} 


// ---
// Bernstein "divide by" operator:
//
// The 'result' vector contains both the quotient
// and remainder vectors in the following format:
//  [result] = [quotient|remainder]
//
// One can use the functions quo() or rem() defined
// later to further extract the quotient or the 
// remainder parts, respectively. 

Bernstein Bernstein::operator/(const Bernstein &uu) const
{
  int ii, jj, jmin, jmax, kk, mm, nn, chgTo;
  mm = dgr;  nn = uu.dgr;
  double **HA, BINOM, pivot, buffer;
  Bernstein result, e1, e2, e12, e3, e123;

  // Set up an m+1 by m+2 matrix:
  try {  
    HA = new double *[mm+1];
    for ( ii = 0; ii <= mm; ii++ )
      *(HA+ii) = new double[mm+2]; 
  }
  catch ( bad_alloc exception ) {  
    cout << "In function divide(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  result.dgr = mm;
  try {  result.cf = new double[result.dgr + 1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function divide(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( ii = 0; ii <= result.dgr; ii++ ) { 
    result.cf[ii] = 0.0;
    for ( jj = 0; jj <= mm + 1; jj++ )  
      HA[ii][jj] = 0.0;
  }

  // Construct the big [H|A] matrix:
  for ( kk = 0; kk <= mm; kk++ ) {
    jmin = kk - nn;  
    if ( jmin < 0 )    jmin = 0;
    jmax = kk;      
    if ( jmax > mm-nn )  jmax = mm - nn;

    for ( jj = jmin; jj <= jmax; jj++ ) {  
      e1 = factorization(mm-nn, jj);
      e2 = factorization(nn, kk-jj);
      e3 = factorization(mm, kk);  
      e12 = binomult(e1, e2);      
      e123 = binodiv(e12, e3);    
      BINOM = binom(e123);
      HA[kk][jj] = BINOM * uu.cf[kk-jj];
    }

    jmin = kk - mm + nn -1;  
    if ( jmin < 0 )    jmin = 0;
    jmax = kk;             
    if ( jmax > nn-1 )  jmax = nn - 1;

    for ( jj = jmin; jj <= jmax; jj++ ) {
      e1 = factorization(mm-nn+1, kk-jj);  
      e2 = factorization(nn-1, jj);    
      e3 = factorization(mm, kk);      
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);
      BINOM = binom(e123);        
      HA[kk][mm-nn+1+jj] = BINOM;
    }
    HA[kk][mm+1] = cf[kk]; 
  }

  // Gaussian elimination with partial pivoting:
  for ( ii = 0; ii < mm; ii++ ) {
    // Interchanging rows if necessary:
    pivot = HA[ii][ii]; 
    chgTo = 0;
    for ( jj = ii + 1; jj <= mm; jj++ ) { 
      if ( fabs(pivot) < fabs(HA[jj][ii]) )  { 
        pivot = HA[jj][ii]; 
        chgTo = jj; 
      }
    }
    if ( pivot != HA[ii][ii] ) {
      // Interchange two rows:
      for ( kk = 0; kk <= mm + 1; kk++ ) { 
        buffer = HA[ii][kk]; 
        HA[ii][kk] = HA[chgTo][kk]; 
        HA[chgTo][kk] = buffer;  
      }
    }
    // LU factorization:
    for ( jj = ii + 1; jj <= mm; jj++ ) {
      // Compute and store the multipliers:
      HA[jj][ii] = HA[jj][ii] / HA[ii][ii];    

      for ( kk = ii + 1; kk <= mm + 1; kk++ )
        HA[jj][kk] = HA[jj][kk] - HA[jj][ii] * HA[ii][kk];
    } 
  }

  // back substitution:
  result.cf[mm] = HA[mm][mm+1] / HA[mm][mm];
  for ( ii = mm - 1; ii >= 0; ii-- ) {
    for ( jj = ii + 1; jj <= mm; jj++ )
      result.cf[ii] = result.cf[ii] - HA[ii][jj] * result.cf[jj]; 

    result.cf[ii] = result.cf[ii] + HA[ii][mm+1];
    result.cf[ii] = result.cf[ii] / HA[ii][ii];
  }

  for ( ii = 0; ii <= mm; ii++ )
    delete [] *(HA+ii);
  delete [] HA;
  return result;
}

// ---
// Extract quotient after division:
Bernstein quo(const Bernstein &uu, int degree)
{
  int ii;
  Bernstein result;
  
  result.dgr = degree;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( ii = 0; ii <= result.dgr; ii++ )
    result.cf[ii] = uu.cf[ii];

  return result;
}

// ---
// Extract remainder after division:
Bernstein rem(const Bernstein &uu, int degree)
{
  int ii;
  Bernstein result;

  result.dgr = degree;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function rem(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }


  for ( ii = 0; ii <= result.dgr; ii++ )
    result.cf[ii] = uu.cf[uu.dgr-result.dgr+ii];

  return result;
}

// ---
// Bernstein "composition" operator:
Bernstein Bernstein::operator<<(const Bernstein &uu) const
{
  int ss, ii, jj, kk, ll, kmin, kmax, mm, nn;
  double ***hh, fracBINOM;
  Bernstein result, e1, e2, e12, e3, e123, CF;

  mm = dgr;  nn = uu.dgr;
  result.dgr = mm * nn;

  try {
    result.cf = new double[result.dgr + 1];
    hh = new double **[mm+1];
    for( ii = 0; ii <= mm; ii++ )
      *(hh + ii) = new double *[mm+1];
    for( ii = 0; ii <= mm; ii++ )
      for( jj = 0; jj <= mm; jj++ )
        *(*(hh + ii) + jj) = new double [mm*nn+1];
  }
  catch ( bad_alloc exception ) {  
    cout << "In function composition(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( ii = 0; ii <= mm; ii++ )
    hh[0][ii][0] = cf[ii];

  for ( ss = 1; ss <= mm; ss++ ) {
    for ( ii = 0; ii <= mm - ss; ii++ ) {
      for ( jj = 0; jj <= nn * ss; jj++ ) {
        kmin = jj - nn;      
        if ( kmin < 0 )   kmin = 0;
        kmax = nn * ss - nn;  
        if ( kmax > jj )   kmax = jj;
        hh[ss][ii][jj] = 0.0;

        // Compute the common factor:
        e1 = factorization(nn*ss-nn, kmin);
        e2 = factorization(nn, jj-kmin);
        e3 = factorization(nn*ss, jj); 
        e12 = binomult(e1, e2);
        e123 = binodiv(e12, e3);
        CF = e123;

        for ( kk = kmin + 1; kk <= kmax; kk++ ) {
          e1 = factorization(nn*ss-nn, kk);
          e2 = factorization(nn, jj-kk);
          e3 = factorization(nn*ss, jj); 
          e12 = binomult(e1, e2);
          e123 = binodiv(e12, e3);
     
          for ( ll = 1; ll <= e123.dgr; ll++ ) {
            if ( e123.cf[ll] != CF.cf[ll] )
              CF.cf[ll] = 0;
	  }
	}

	// Bring out the common factor:
        for ( kk = kmin; kk <= kmax; kk++ ) {
          e1 = factorization(nn*ss-nn, kk);
          e2 = factorization(nn, jj-kk);
          e3 = factorization(nn*ss, jj); 
          e12 = binomult(e1, e2);
          e123 = binodiv(e12, e3);

          for ( ll = 1; ll <= e123.dgr; ll++ ) {
            if ( e123.cf[ll] == CF.cf[ll] )
              e123.cf[ll] = 0;
	  }

        fracBINOM = binom(e123);
        hh[ss][ii][jj] += fracBINOM * 
		((1-uu.cf[jj-kk])*hh[ss-1][ii][kk] + uu.cf[jj-kk]*hh[ss-1][ii+1][kk]);   
	}

        hh[ss][ii][jj] = hh[ss][ii][jj] * binom(CF);
      }
    }
  }

  for ( jj = 0; jj <= result.dgr; jj++ )
    result.cf[jj] = hh[mm][0][jj]; 

  for ( ii = 0; ii <= mm; ii++ )
    for ( jj = 0; jj <= mm; jj++ )
      delete [] *(*(hh+ii)+jj);

  for( ii = 0; ii <= mm; ii++ )
    delete [] *(hh+ii);

  delete [] hh;
  return result;
}

// ---
// Bernstein "power" operator:
//
// It only accepts an integer for the power
Bernstein Bernstein::operator^(int power)
{
  int ii;
  Bernstein result;

  result.dgr = 0;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function power(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }
  result.cf[0] = 1.0;

  if ( power == 0 )
    return result;

  for ( ii = 1; ii <= power; ii++ )
    result = result * (*this);

  return result;
}

// ---
// Bernstein "pre-multiply" operator:
Bernstein operator*(double factor, const Bernstein &uu) 
{
  int kk;
  Bernstein result;

  result.dgr = uu.dgr;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function pre-multiply(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( kk = 0; kk <= result.dgr; kk++ )
    result.cf[kk] = factor * uu.cf[kk];    
 
 return result;
}

// ---   
// Bernstein "post-multiply" operator:
Bernstein operator*(const Bernstein &uu, double factor) 
{
  int kk;
  Bernstein result;

  result.dgr = uu.dgr;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( kk = 0; kk <= result.dgr; kk++ )
    result.cf[kk] = factor * uu.cf[kk];    

  return result;
}  

// ---
// Bernstein "differentiate" function:
Bernstein diff(const Bernstein &uu)
{
  int kk; 
  Bernstein result;

  result.dgr = uu.dgr - 1;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( kk = 0; kk <= result.dgr; kk++ )
    result.cf[kk] = (double)uu.dgr * (uu.cf[kk+1] - uu.cf[kk]);

  return result;
}  

// ---
// Bernstein "integrate" (indefinite) function:
Bernstein integrate(const Bernstein &uu)
{
  int kk, jj; 
  double recip = 1.0 / (double)(uu.dgr+1);
  Bernstein result;

  result.dgr = uu.dgr + 1;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( kk = 0; kk <= result.dgr; kk++ )
    result.cf[kk] = 0.0;
  for ( kk = 1; kk <= result.dgr; kk++ )
    for ( jj = 0; jj <= kk-1; jj++ )
      result.cf[kk] += recip * uu.cf[jj];

  return result;
}

// ---
// Bernstein "integrate" (definite) function:
double integral(const Bernstein &uu)
{
  int kk, jj; 
  double recip = 1.0 / (double)(uu.dgr+1);
  Bernstein result;

  result.dgr = uu.dgr + 1;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    return 0.0;
  }

  for ( kk = 0; kk <= result.dgr; kk++ )
    result.cf[kk] = 0.0;
  for ( kk = 1; kk <= result.dgr; kk++ )
    for ( jj = 0; jj <= kk - 1; jj++ )
      result.cf[kk] += recip * uu.cf[jj];

  return result.cf[result.dgr];
}

// ---
// Function for finding the L2 norm:
double Norm(const Bernstein &uu)
{
  int ii, jj, nn;
  double BINOM, Norm = 0.0;
  Bernstein e1, e2, e12, e3, e123;
  nn = uu.dgr;

  // There's no need to bring out the common factors because it's 1.
  for ( ii = 0; ii <= nn; ii++ ) {
    for ( jj = 0; jj <= nn; jj++ ) {
      e1 = factorization(nn, ii);
      e2 = factorization(nn, jj);
      e3 = factorization(nn+nn, ii+jj); 
      e12 = binomult(e1, e2);
      e123 = binodiv(e12, e3);
      BINOM = binom(e123);
      Norm += BINOM * uu.cf[ii] * uu.cf[jj];
    }
  }

  Norm = Norm / (2.0 * (double)nn + 1.0);
  Norm = sqrt(Norm);
  return Norm;
}

// ---
// Function for normalizing a polynomial in Bernstein form:
Bernstein Normalize(const Bernstein &uu)
{
  int ii, nn;
  double factor;  
  Bernstein result;
  nn = uu.dgr;

  result.dgr = nn;
  try {  result.cf = new double[result.dgr+1];  }
  catch ( bad_alloc exception ) {  
    cout << "In function quo(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  factor = 1.0 / Norm(uu);
  for ( ii = 0; ii <= uu.dgr; ii++ )
    result.cf[ii] = factor * uu.cf[ii];

  return result;
}

// ---
// Bernstein "greatest common divisor" function:
Bernstein GCD(const Bernstein &uu, const Bernstein &vv, double epsilon) 
{
  int nn;
  double r1Norm, r2Norm;
  Bernstein ff, gg, phi0, phi1, phi2, result; 
 
  // Initialization:
  if ( uu.dgr >= vv.dgr ) { 
    nn = vv.dgr;
    ff = Normalize(uu);     
    gg = Normalize(vv);
  }
  else { 
    nn = uu.dgr; 
    ff = Normalize(vv);     
    gg = Normalize(uu);
  }

  phi0 = ff;  phi1 = gg;

  if ( nn == 0 )
    result = phi1;

  while ( nn > 0 )
  {
    // Divide, obtain the remainder and store it in phi2:
    phi2 = rem( phi0 / phi1, nn - 1 );  

    // Compute the sizes of the remainders, when
    // the (normalized) original polynomials are
    // divided by the `gcd candidate' phi1:
    r1Norm = Norm( rem( ff / phi1, phi1.dgr - 1 ) );
    r2Norm = Norm( rem( gg / phi1, phi1.dgr - 1 ) );
    //  printf("divisor degree = %d   r1Norm = %1.16f\n",phi1.dgr,r1Norm);
    //  printf("divisor degree = %d   r2Norm = %1.16f\n",phi1.dgr,r2Norm);
    if ( r1Norm <= epsilon && r2Norm <= epsilon ) {
      // Found gcd at the mth stage:
      break;
    }   

    // If no gcd is found, proceed to the next division:
    phi0 = phi1;
    phi1 = phi2;  
    nn = phi2.dgr;  
  }

  result = phi1;  
  return result;
}

// ---
// Function for sorting the cf[]'s of a Bernstein object, 
// in ascending order.
// This function is used with ROOT() to sort the roots found.
Bernstein sort(const Bernstein &uu) 
{
  int ii, jj, nn, count, minindex;
  double min, *temp, *newtemp;
  Bernstein result;
 
  nn = uu.dgr;
  result.dgr = nn;
  if ( result.dgr == -1 )
    return result;

  try {  
    // Beware: n elements instead of n+1 !!
    result.cf = new double[result.dgr];  
    temp = new double[result.dgr];
  }
  catch ( bad_alloc exception ) {  
    cout << "In function sort(): " << exception.what() << endl;
    result.dgr = -1;
    return result;
  }

  for ( ii = 0; ii < result.dgr; ii++ )
    temp[ii] = uu.cf[ii];

  for ( ii = result.dgr; ii >= 2; ii-- ) { 
    min = temp[ii-1];  
    minindex = ii - 1;
    for ( jj = ii - 2; jj >= 0; jj-- ) {
      if ( temp[jj] <= min ) {  
        min = temp[jj];  
        minindex = jj;  
      }
    }

    result.cf[result.dgr-ii] = min;

    try {  newtemp = new double[ii-1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sort(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    count = 0;
    for ( jj = 0; jj < ii; jj++ ) {    
      if ( jj != minindex ) {  
        newtemp[count] = temp[jj]; 
        count += 1;  
      }
    }

    delete [] temp;

    try {  temp = new double[ii-1];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sort(): " << exception.what() << endl;
      result.dgr = -1;
      return result;
    }

    for ( jj = 0; jj < ii - 1; jj++ )
      temp[jj] = newtemp[jj];

    delete [] newtemp;
  }

  result.cf[nn-1] = temp[0];
  return result;
}


// ---
// Bernstein root solver
//
// epsilon is the tolerance for gcd()
//
// delta is the tolerance defined as: 
// Whenever |F(u)| <= delta, u is identified as a root.
//
// eta is the tolerance defined as:
// Whenever the width of a sub-interval b - a < eta, we
// don't subdivide it anymore --- one can think of it 
// as the resolution of subdivisions.
Bernstein ROOT(const Bernstein &uu, double delta, double eta, double epsilon)
{
  int ii, mm, count;	
  Bernstein uprime, DD, roots;

  FILE *ROOTlog;
  ROOTlog = fopen("ROOT.log","w");

  uprime = diff(uu);
  DD = GCD(uu,uprime,epsilon);

  if ( DD.dgr == 0 ) 
    // There are no multiple roots: invoke the simple root solver
    roots = sROOT(uu, delta, eta);  

  else
    // There are multiple roots: invoke the multiple root solver
    roots = mROOT(uu, delta, eta, epsilon);  

  mm = 1;  ii = 0;  count = 0;

  while( (ii+1) < roots.dgr ) {
    if ( roots.cf[ii+1] != roots.cf[ii] ) {  
      count++;
      fprintf(ROOTlog,
          "Root %d (multiplicity %d): %1.6f,\n",count,mm,roots.cf[ii]);
      mm = 1;
    }

    else
      mm++;

    ii++; 
  }

  count++;
  fprintf(ROOTlog,"Root %d (multiplicity %d): %1.6f,\n",count,mm,roots.cf[ii]);

  fprintf(ROOTlog,"\nTotal of %d distinct roots are found.\n",count);

  return roots;
}

// ---
// Bernstein root solver for simple roots: 
// This function returns a Bernstein instance for which
// dgr = number of roots found; and
// cf[0],...,cf[dgr-1] record these roots    
Bernstein sROOT(const Bernstein &uu, double delta, double eta)
{
  int ii, jj, kk, nn, last_GO, GO, OUT, abandon, ONLY_ONE, ALL_ONE;
  int *signchg, *GO_index, *OUT_index, converge;
  double step, base, min2, max2, min3, max3; 
  double root, error, iterate;
  Bernstein uprime, yprime, ydbprime, roots, temp;
  Bernstein *x_stack, *y_stack, *last_x, *last_y, *ONE_x;
  FILE *sROOTlog;

  sROOTlog = fopen("sROOT.log", "w");

  nn = uu.dgr;
  uprime = diff(uu);
  // Prepare root storage array:
  roots.dgr = 0;
  try {  roots.cf = new double[nn];  }  
  // Beware: roots.cf[n] is NOT used !!
  catch ( bad_alloc exception ) {  
    cout << "In function sROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }

  for ( ii = 0; ii < nn; ii++ )
    roots.cf[ii] = -1;  

  // Initialize:
  //
  // The variable `GO' keeps track of the number of 
  // sub-sections that the algorithm is processing.
  // The variable `OUT' keeps track of the number of 
  // sub-sections that pass the tests, and are thus 
  // discarded from `GO'. 
  ALL_ONE = 0;
  GO = 1;  OUT = 0;
  ONLY_ONE = 0;

  try {  
    last_x = new Bernstein[GO];
    last_y = new Bernstein[GO];
    // Assign the starting Bernstein data (x,y):
    last_x[0].dgr = nn;
    last_x[0].cf = new double[nn+1];
    // Prepare the array `ONE_x' for storing the sub-intervals 
    // that have exactly one coefficient sign change:
    ONE_x = new Bernstein[2*nn];
  }
  catch ( bad_alloc exception ) {  
    cout << "In function sROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }
 
  step = 1.0 / (double)nn; 
  for ( ii = 0; ii <= nn; ii++ )
    last_x[0].cf[ii] = 0.0 + (double)ii * step;

  // Normalize the original polynomial:
  last_y[0] = Normalize(uu);

  // When there are still sub-sections that have more than 1 
  // x-axis intersections, continue this while loop:
  while(ALL_ONE == 0)
  {
    last_GO = GO;  // number of sub-sections conveyed from last step 
    try {
      GO_index = new int[2*last_GO];
      signchg = new int[2*last_GO];
    }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    // Examine subdivided sections. Mark the sub-sections
    // that "pass" the tests with "GO":

    // Test 1. Check sign changes:
    GO = 0;  // Reset the "GO" counter.
    for ( ii = 0; ii < last_GO; ii++ ) 
    {
      // Examine the left half:
      temp = subLEFT(last_y[ii],0.5);  
      signchg[2*ii] = 0;
      for ( jj = 0; jj < nn; jj++ )
      {
        if ((temp.cf[jj] * temp.cf[jj+1]) <= 0)
          signchg[2*ii] += 1;
        else if ( jj == 0 && (temp.cf[0] < delta && temp.cf[0] > -1.0*delta) )
          signchg[2*ii] += 1;
        else if ( jj == nn-1 && (temp.cf[nn] < delta && temp.cf[nn] > -1.0*delta) )
          signchg[2*ii] += 1;
      }
      if (signchg[2*ii] > 0)
      {
        GO = GO + 1;  // Convey the sub-sections that
                      // have at least 1 x-axis intersections.
        GO_index[GO-1] = 2*ii;
      }
   
      // Examine the right half:
      temp = subRIGHT(last_y[ii],0.5);
      signchg[2*ii+1] = 0;
      for ( jj = 0; jj < nn; jj++ )
      {
        if (temp.cf[jj] * temp.cf[jj+1] <= 0)
          signchg[2*ii+1] += 1;
        else if ( jj == 0 && (temp.cf[0] < delta && temp.cf[0] > -1.0*delta)  )
          signchg[2*ii+1] += 1;
        else if ( jj == nn-1 && (temp.cf[nn] < delta && temp.cf[nn] > -1.0*delta)  )
          signchg[2*ii+1] += 1;
      }
      if (signchg[2*ii+1] > 0)
      {
        GO = GO + 1;  // Convey the sub-sections that
                      // have at least 1 x-axis intersections.
        GO_index[GO-1] = 2*ii+1;
      }
    }

    // Subdivide according to "GO_index":
    // (only keep sub-sections of concern to save memory) 
    try {
      x_stack = new Bernstein[GO];
      y_stack = new Bernstein[GO];
    }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    } 

    for ( ii = 0; ii < GO; ii++ )
    {
      if ((GO_index[ii]+1)/2 == (GO_index[ii])/2)
      // If it's an even number, it's a left half:
      {
        x_stack[ii] = subLEFT(last_x[GO_index[ii]/2],0.5);
        y_stack[ii] = subLEFT(last_y[GO_index[ii]/2],0.5);  
      }
      else
      // Otherwise it's a right half:
      {
        x_stack[ii] = subRIGHT(last_x[GO_index[ii]/2],0.5);
        y_stack[ii] = subRIGHT(last_y[GO_index[ii]/2],0.5);  
      } 
    }

    delete [] last_x;  delete [] last_y;  
    delete [] signchg;  delete [] GO_index;
  
    // Beyond this point all of the sub-sections conveyed to 
    // the next steps have at least 1 intersection with the 
    // x-axis!!

    fprintf(sROOTlog,
        "Sub-intervals found that have at least 1 x-axis intersections:\n");
    for ( ii = 0; ii < GO; ii++ )
      fprintf(sROOTlog,"[%1.7f, %1.7f]\n",x_stack[ii].cf[0],x_stack[ii].cf[nn]);

    // Test 1. Check if F(a)F(b) < 0
    // Test 2. Check if F'(x) != 0 for x in [a,b]
    // Test 3. Check if F''(x) either >= 0 or <= 0
    OUT = 0; 
    try {  OUT_index = new int[GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    for ( ii = 0; ii < GO; ii++ )
    {
      yprime = diff(y_stack[ii]);  
      ydbprime = diff(yprime);

      max2 = yprime.cf[0];
      for ( kk = 1; kk <= nn - 1; kk++ )
      {  if (yprime.cf[kk] >= max2)   max2 = yprime.cf[kk];  }

      min2 = yprime.cf[0];
      for ( kk = 1; kk <= nn - 1; kk++ )
      {  if (yprime.cf[kk] <= min2)   min2 = yprime.cf[kk];  }

      max3 = ydbprime.cf[0];
      for (kk = 1; kk <= nn - 2; kk++ )
      {  if (ydbprime.cf[kk] >= max3)   max3 = ydbprime.cf[kk];  }

      min3 = ydbprime.cf[0];
      for ( kk = 1; kk <= nn - 2; kk++ )
      {  if (ydbprime.cf[kk] <= min3)   min3 = ydbprime.cf[kk];  }

      if ((y_stack[ii].cf[0]*y_stack[ii].cf[nn] <= delta) &&
		  ((max2*min2) > 0) && ((max3*min3) >= 0))
      {
        // Pass the tests!
        ONLY_ONE = ONLY_ONE + 1;
	// Record the sub-interval that has only 1 x-axis 
        // intersection:
        ONE_x[ONLY_ONE-1] = x_stack[ii];
        OUT = OUT + 1;
        OUT_index[OUT-1] = ii;
      }
    }

    if (OUT != GO)
    { 
      // When there are still sub-sections that have more 
      // than 1 x-axis intersections, continue this while loop:
      ALL_ONE = 0;
      try {
        last_x = new Bernstein[GO-OUT]; 
        last_y = new Bernstein[GO-OUT];
      }
      catch ( bad_alloc exception ) {  
        cout << "In function sROOT(): " << exception.what() << endl;
        roots.dgr = -1;
        return roots;
      }

      last_GO = GO;  GO = 0;

      // Abandon the sub-interval that has only 1 x-axis 
      // intersection.
      // Store the "non-abandoned" intervals:
      for ( ii = 0; ii < last_GO; ii++ ) 
      { 
        abandon = 0;
        for ( jj = 0; jj < OUT; jj++ ) 
        {
          if ( ii == OUT_index[jj] )
          {  abandon = 1;  break;  }  
	}                                 

        if ( abandon == 0 )  
        {
          GO = GO + 1;                    
          last_x[GO-1] = x_stack[ii];
          last_y[GO-1] = y_stack[ii];
        }
      }
    } 

    else
      ALL_ONE = 1;

    delete [] x_stack;  delete [] y_stack;  
    delete [] OUT_index;
  }

  // Beyond this point all of the sub-sections have only
  // 1 intersection with the x-axis:

  fprintf(sROOTlog,"The conveyed sub-intervals after the first stage:\n");
  for ( ii = 0; ii < ONLY_ONE; ii++ )
  fprintf(sROOTlog,"[%1.7f, %1.7f]\n",ONE_x[ii].cf[0],ONE_x[ii].cf[nn]);

  // Prepare for the next while loop:
  GO = ONLY_ONE;  converge = 0;

  try {  x_stack = new Bernstein[GO];  }   
  catch ( bad_alloc exception ) {  
    cout << "In function sROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }

  for ( ii = 0; ii < GO; ii++ )
    x_stack[ii] = ONE_x[ii];

  delete [] ONE_x;    

  while(converge == 0)
  {       
    converge = 1;
    // Test 4. Find endpoint roots and abandon the 
    // corresponding sub-intervals:
    OUT = 0;                   
    try {  OUT_index = new int[GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    fprintf(sROOTlog,"\nEnd-point roots:\n");

    for ( ii = 0; ii < GO; ii++ )
    {  
      if ( fabs(EVAL(uu,x_stack[ii].cf[0])) <= delta )    
      {
        abandon = 0;	 
	for ( jj = 0; jj < roots.dgr; jj++ ) {
          if ( fabs(roots.cf[jj] - x_stack[ii].cf[0]) <= 1e-15 )
            abandon = 1;
	}
	if ( abandon == 0 ) {
          roots.dgr = roots.dgr + 1;      
          roots.cf[roots.dgr - 1] = x_stack[ii].cf[0];
		  fprintf(sROOTlog,"One endpoint root found: %1.16f\n",x_stack[ii].cf[0]);
        }
        OUT = OUT + 1;
        OUT_index[OUT-1] = ii;    
      }

      else if ( fabs(EVAL(uu,x_stack[ii].cf[nn])) <= delta )    
      {
        abandon = 0;	 
	for ( jj = 0; jj < roots.dgr; jj++ ) {
          if ( fabs(roots.cf[jj] - x_stack[ii].cf[nn]) <= 1e-15 )
            abandon = 1;
	}
	if ( abandon == 0 ) {
          roots.dgr = roots.dgr + 1;      
          roots.cf[roots.dgr - 1] = x_stack[ii].cf[nn];
		  fprintf(sROOTlog,"One endpoint root found: %1.16f\n",x_stack[ii].cf[nn]);
	}
        OUT = OUT + 1;
        OUT_index[OUT-1] = ii;		   
      }
    }    

    try {  last_x = new Bernstein[GO-OUT];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    last_GO = GO;
    GO = 0;
    // Store the "non-abandoned" sub-intervals:
    for ( ii = 0; ii < last_GO; ii++ ) 
    { 
      abandon = 0;
      for ( jj = 0; jj < OUT; jj++ ) 
      {
        if ( ii == OUT_index[jj] )
        {  abandon = 1;  break;  }  
      }                                

      if (abandon == 0)  
      {
        GO = GO + 1;                    
        last_x[GO-1] = x_stack[ii];
      }
    }

    delete [] x_stack;  
    delete [] OUT_index;

    if (GO == 0)
      // All roots are found!
      return roots;

    // Test 5. The tangent at the endpoint of the sub-interval 
    //    that has a smaller absolute slope, must intersect 
    //    the x-axis at a point within the sub-interval:
    OUT = 0;                   
    try {  OUT_index = new int[GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    fprintf(sROOTlog,"\nSub-intervals that fail test 5:\n");

    for ( ii = 0; ii < GO; ii++ )
    {              
      if (fabs(EVAL(uprime,last_x[ii].cf[nn])) <= 
          fabs(EVAL(uprime,last_x[ii].cf[0])))
        base = EVAL(uu,last_x[ii].cf[nn]) / EVAL(uprime,last_x[ii].cf[nn]);
      else
        base = EVAL(uu,last_x[ii].cf[0]) / EVAL(uprime,last_x[ii].cf[0]);

      if ((fabs(base) > (last_x[ii].cf[nn] - last_x[ii].cf[0])) &&
          (fabs(last_x[ii].cf[nn] - last_x[ii].cf[0]) >= eta))  
      // When the sub-interval [a,b] is too small, b - a 
      // entails catastrophic cancellations!!

      {  // It won't converge
	 converge = 0;  
	 fprintf(sROOTlog,"[%1.7f, %1.7f]:\n",last_x[ii].cf[0],last_x[ii].cf[nn]);
         fprintf(sROOTlog,"|base| = %1.16f\n",fabs(base)); 
      }

      else
      {    
        // Pass the test, do Newton's iterations:
        // It should converge from arbitrary start point:
        root = last_x[ii].cf[nn/2];
        iterate = 50;
        for ( jj = 0; jj <= iterate; jj++ )
        {
          error = EVAL(uu,root) / EVAL(uprime,root);
          root = root - error;
          if ( fabs(EVAL(uu,root)) <= delta )
          {
            roots.dgr = roots.dgr + 1; 
            roots.cf[roots.dgr - 1] = root;
	    jj = 100;
            break;
          }
        }
      
        if (jj != 100)
        {
          // If it doesn't converge to meet the root tolerance:
          fprintf(sROOTlog,
             "Warning: This value %1.16f does not yield convergence!\n",root);  
	  fprintf(sROOTlog,"value = %1.16f\n",fabs(EVAL(uu,root)));
	  // But record it anyways, since there's not much
	  // we can do:
          roots.dgr = roots.dgr + 1; 
          roots.cf[roots.dgr - 1] = root;  
	}

        OUT = OUT + 1;
        OUT_index[OUT-1] = ii;
      }
    }

    // Abandon the unwanted intervals:
    try {  x_stack = new Bernstein[GO-OUT];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    last_GO = GO;
    GO = 0;
    // Store the "non-abandoned" intervals:
    for ( ii = 0; ii < last_GO; ii++ ) 
    { 
      abandon = 0;
      for ( jj = 0; jj < OUT; jj++ ) 
      {
        if ( ii == OUT_index[jj] )
        {  abandon = 1;  break;  }  
      }                               

      if ( abandon == 0 )  
      {
        GO = GO + 1;                    
        x_stack[GO-1] = last_x[ii];
      }
    }

    delete [] last_x;  
    delete [] OUT_index;

    if (GO == 0)
    // All roots are found!
      return roots;
   
    // Examine subdivided sections. Mark the intervals
    // that "pass" the tests with "GO":
  
    // Check sign changes in the interval:
    last_GO = GO;  
    try {  GO_index = new int[2*last_GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }

    GO = 0;  // Reset the "GO" counter.
    for ( ii = 0; ii < last_GO; ii++ ) 
    {
      // Examine the left half:
      temp = subLEFT(x_stack[ii],0.5);
    
      if ( EVAL(uu,temp.cf[0]) * EVAL(uu,temp.cf[nn]) < 0 )
      {   
        GO = GO + 1;  // Abandon the sub-intervals that
                      // don't have exactly one x-axis intersections.
        GO_index[GO-1] = 2*ii;
      }
   
      // Examine the right half:
      temp = subRIGHT(x_stack[ii],0.5);

      if ( EVAL(uu,temp.cf[0]) * EVAL(uu,temp.cf[nn]) < 0 )
      {   
        GO = GO + 1;  // Abandon the sub-intervals that
                      // don't have exactly one x-axis intersections.
        GO_index[GO-1] = 2*ii+1;
      }
    }
     
    // Subdivide according to "GO_index":
    try {  last_x = new Bernstein[GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }
   
    for ( ii = 0; ii < GO; ii++ )
    {
      if ((GO_index[ii]+1)/2 == (GO_index[ii])/2)
        // If it's an even number, it's a left half:
        last_x[ii] = subLEFT(x_stack[GO_index[ii]/2],0.5);
      else
        // Otherwise it's the right half:
        last_x[ii] = subRIGHT(x_stack[GO_index[ii]/2],0.5);
    }

    delete [] x_stack;    delete [] GO_index;

    try {  x_stack = new Bernstein[GO];  }
    catch ( bad_alloc exception ) {  
      cout << "In function sROOT(): " << exception.what() << endl;
      roots.dgr = -1;
      return roots;
    }
   
    for ( ii = 0; ii < GO; ii++ )
      x_stack[ii] = last_x[ii];    

    delete [] last_x;  
  } 
  return roots;
}

// ---
// Bernstein root solver for multiple roots:
// epsilon is the tolerance for gcd()
Bernstein mROOT(const Bernstein &uu, double delta, double eta, double epsilon)
{
  int nn, kk, mm, ii, jj;
  Bernstein up, DD, Dp, *Dk, *ff, *XX;
  Bernstein roots, temp;

  nn = uu.dgr;
  // Prepare root storage array:
  try {  roots.cf = new double[nn];  }  
  // Beware: roots.cf[n] is NOT used !!
  catch ( bad_alloc exception ) {  
    cout << "In function mROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }

  for ( ii = 0; ii < nn; ii++ )
    roots.cf[ii] = -1;

  // Initialize:
  up = diff(uu);

  try {  Dk = new Bernstein[uu.dgr];  }
  catch ( bad_alloc exception ) {  
    cout << "In function mROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }

  Dk[0] = GCD(uu,up,epsilon);

  // Compute D's and determine m:
  kk = 1;
  while (1) {
    Dp = diff(Dk[kk-1]); 
    Dk[kk] = GCD(Dk[kk-1],Dp,epsilon);

    if (Dk[kk].dgr == 0) {
      mm = kk + 1;  // no roots of multiplicity higher than m
      break;
    }
    kk = kk + 1; 
  }
 
  // Compute f's:
  try {  
    ff = new Bernstein[mm + 1];  
    XX = new Bernstein[mm];
  }
  catch ( bad_alloc exception ) {  
    cout << "In function mROOT(): " << exception.what() << endl;
    roots.dgr = -1;
    return roots;
  }

  ff[0] = uu;
  ff[1] = quo ( ff[0] / Dk[0], ff[0].dgr - Dk[0].dgr );
  for ( ii = 2; ii <= mm; ii++ )
    ff[ii] = quo ( Dk[ii-2] / Dk[ii-1], Dk[ii-2].dgr - Dk[ii-1].dgr );

  delete [] Dk;

  // Compute X's, solve simple roots and record them
  // according to respective multiplicities:
  for ( ii = 0; ii < mm - 1; ii++ ) {
    XX[ii] = quo ( ff[ii+1] / ff[ii+2], ff[ii+1].dgr - ff[ii+2].dgr );
    temp = sROOT(XX[ii], delta, eta);
    for ( jj = 0; jj < (temp.dgr*(ii+1)); jj++ ) {
      roots.dgr = roots.dgr + 1;
      roots.cf[roots.dgr-1] = temp.cf[jj/(ii+1)];
    }
  }

  XX[mm-1] = ff[mm];
  temp = sROOT(XX[mm-1], delta, eta);
  for ( jj = 0; jj < (temp.dgr*mm); jj++ )
  {
    roots.dgr = roots.dgr + 1;
    roots.cf[roots.dgr-1] = temp.cf[jj/mm];
  }

  delete [] ff;  delete [] XX;
  return roots;
}

















































