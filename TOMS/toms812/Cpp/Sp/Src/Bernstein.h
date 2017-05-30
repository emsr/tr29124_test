#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <new>
using std::bad_alloc;

// Bernstein.h:
// Declaration of a class for computations with polynomials 
// in Bernstein form.
// 
// by Evan Yi-Feng Tsai, 2001

class Bernstein
{
 public:   
   double *cf;  // Array for storing Bernstein coefficients
                // (It also stores exponents of the factorization
                //  of binomial coefficients).
   int dgr;     // Degree of the Bernstein polynomial
                // (Also the number of elements of the 
                //  exponent array).
   
   // constructors and destructor: 
   Bernstein();                           // default constructor 
   Bernstein(double *coeff, int degree);  // constructor with coefficient 
                                          // and degree inputs
   Bernstein(const Bernstein &u);         // copy constructor   
   ~Bernstein();                          // destructor

   // Binomial functions and utilities:
   friend Bernstein getprimes(int n, int pr_num, int mode);
   friend Bernstein factorization(int n, int k);
   friend double binom(int n, int k);   
   friend double binom(const Bernstein &e);
   friend Bernstein binomult(const Bernstein &e1, const Bernstein &e2);
   friend Bernstein binodiv(const Bernstein &e1, const Bernstein &e2);

   // Degree elevation:
   friend Bernstein DE(const Bernstein &u, int r);

   // Evaluation and subdivision:
   friend double EVAL(const Bernstein &u, double x);
   friend Bernstein subLEFT(const Bernstein &u, double x);
   friend Bernstein subRIGHT(const Bernstein &u, double x);

   // Differentiation and integration:
   friend Bernstein diff(const Bernstein &u);
   friend Bernstein integrate(const Bernstein &u);
   friend double integral(const Bernstein &u);

   // Normalization:
   friend double Norm(const Bernstein &u);
   friend Bernstein Normalize(const Bernstein &u);

   // Operator overloading:
   const Bernstein &operator=(const Bernstein &u);
   Bernstein operator+(const Bernstein &u) const;
   Bernstein operator-(const Bernstein &u) const;    
   Bernstein operator*(const Bernstein &u) const;
   friend Bernstein operator*(double factor, const Bernstein &u);
   friend Bernstein operator*(const Bernstein &u, double factor);
   Bernstein operator^(int power);
   Bernstein operator/(const Bernstein &u) const;
   friend Bernstein quo(const Bernstein &u, int degree);
   friend Bernstein rem(const Bernstein &u, int degree);
   Bernstein operator<<(const Bernstein &u) const;

   // Greatest common divisors:
   friend Bernstein GCD(const Bernstein &u, const Bernstein &v, double epsilon);

   // Bernstein root solver:
   friend Bernstein sort(const Bernstein &u);
   friend Bernstein ROOT(const Bernstein &u, double delta, double eta, double epsilon);
   friend Bernstein sROOT(const Bernstein &u, double delta, double eta);
   friend Bernstein mROOT(const Bernstein &u, double delta, double eta, double epsilon);

};















