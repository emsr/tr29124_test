# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "bernstein_polynomial.hpp"

//****************************************************************************80

double *bernstein_matrix ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_MATRIX returns the Bernstein matrix.
//
//  Discussion:
//
//    The Bernstein matrix of order N is an NxN matrix A which can be used to
//    transform a vector of power basis coefficients C representing a polynomial 
//    P(X) to a corresponding Bernstein basis coefficient vector B:
//
//      B = A * C
//
//    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
//    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
//    For N = 5, the matrix has the form:
//
//      1 -4   6  -4  1
//      0  4 -12  12 -4
//      0  0   6 -12  6
//      0  0   0   4 -4
//      0  0   0   0  1
//
//    and the numbers in each column represent the coefficients in the power
//    series expansion of a Bernstein polynomial, so that 
//
//      B(5,4) = - 4 x^4 + 12 x^3 - 12 x^2 + 4 x
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_MATRIX[N*N], the Bernstein matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];
 
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*n] = r8_mop ( j - i ) * r8_choose ( n - 1 - i, j - i ) 
        * r8_choose ( n - 1, i );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double bernstein_matrix_determinant ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_MATRIX_DETERMINANT returns the determinant of the BERNSTEIN matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_MATRIX_DETERMINANT, the determinant.
//
{
  int i;
  double value;

  value = 1.0;
  for ( i = 0; i < n; i++ )
  {
    value = value * r8_choose ( n - 1, i );
  }

  return value;
}
//****************************************************************************80

double *bernstein_matrix_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
//
//  Discussion:
//
//    The inverse Bernstein matrix of order N is an NxN matrix A which can 
//    be used to transform a vector of Bernstein basis coefficients B
//    representing a polynomial P(X) to a corresponding power basis 
//    coefficient vector C:
//
//      C = A * B
//
//    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
//    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
//    For N = 5, the matrix has the form:
//
//      1   1    1    1   1
//      0  1/4  1/2  3/4  1
//      0   0   1/6  1/2  1
//      0   0    0   1/4  1
//      0   0    0    0   1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_MATRIX_INVERSE[N*N], the inverse Bernstein matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*n] = r8_choose ( j, i ) / r8_choose ( n - 1, i );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *bernstein_poly_01 ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01 evaluates the Bernstein polynomials based in [0,1].
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula is:
//
//      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials 
//    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
//    each of degree N, which form a basis for polynomials on [0,1].
//
//    Input, double X, the evaluation point.
//
//    Output, double BERNSTEIN_POLY[N+1], the values of the N+1 
//    Bernstein polynomials at X.
//
{
  double *bern;
  int i;
  int j;

  bern = new double[n+1];

  if ( n == 0 )
  {
    bern[0] = 1.0;
  }
  else if ( 0 < n )
  {
    bern[0] = 1.0 - x;
    bern[1] = x;
 
    for ( i = 2; i <= n; i++ )
    {
      bern[i] = x * bern[i-1];
      for ( j = i - 1; 1 <= j; j-- )
      {
        bern[j] =         x   * bern[j-1] 
                + ( 1.0 - x ) * bern[j];
      }
      bern[0] = ( 1.0 - x ) * bern[0];
    }
  }
  return bern;
}
//****************************************************************************80

double *bernstein_poly_01_matrix ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01_MATRIX evaluates the Bernstein polynomials based in [0,1].
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula is:
//
//      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the degree of the Bernstein polynomials 
//    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
//    each of degree N, which form a basis for polynomials on [0,1].
//
//    Input, double X[M], the evaluation points.
//
//    Output, double BERNSTEIN_POLY_01_MATRIX[M*(N+1)], the values of the N+1 
//    Bernstein polynomials at the evaluation points.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    if ( n == 0 )
    {
      b[i+0*m] = 1.0;
    }
    else if ( 0 < n )
    {
      b[i+0*m] = 1.0 - x[i];
      b[i+1*m] = x[i];
 
      for ( j = 2; j <= n; j++ )
      {
        b[i+j*m] = x[i] * b[i+(j-1)*m];
        for ( k = j - 1; 1 <= k; k-- )
        {
          b[i+k*m] =         x[i]   * b[i+(k-1)*m] 
                   + ( 1.0 - x[i] ) * b[i+k*m];
        }
        b[i+0*m] = ( 1.0 - x[i] ) * b[i+0*m];
      }
    }
  }
  return b;
}
//****************************************************************************80

void bernstein_poly_01_values ( int &n_data, int &n, int &k, double &x, 
  double &b )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01_VALUES returns some values of the Bernstein polynomials.
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula for the Bernstein polynomials is
//
//      B(N,I)(X) = [N!/(I!(N-I)!)] * (1-X)^(N-I) * X^I
//
//    In Mathematica, the function can be evaluated by:
//
//      Binomial[n,i] * (1-x)^(n-i) * x^i
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2  * X
//    B(3,2)(X) = 3 * (1-X)    * X^2
//    B(3,3)(X) =                X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3  * X
//    B(4,2)(X) = 6 * (1-X)^2  * X^2
//    B(4,3)(X) = 4 * (1-X)    * X^3
//    B(4,4)(X) =                X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the degree of the polynomial.
//
//    Output, int &K, the index of the polynomial.
//
//    Output, double &X, the argument of the polynomial.
//
//    Output, double &B, the value of the polynomial B(N,K)(X).
//
{
# define N_MAX 15

  static double b_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.7500000000000000E+00,
     0.2500000000000000E+00,
     0.5625000000000000E+00,
     0.3750000000000000E+00,
     0.6250000000000000E-01,
     0.4218750000000000E+00,
     0.4218750000000000E+00,
     0.1406250000000000E+00,
     0.1562500000000000E-01,
     0.3164062500000000E+00,
     0.4218750000000000E+00,
     0.2109375000000000E+00,
     0.4687500000000000E-01,
     0.3906250000000000E-02 };

  static int k_vec[N_MAX] = {
    0,
    0, 1,
    0, 1, 2,
    0, 1, 2, 3,
    0, 1, 2, 3, 4 };

  static int n_vec[N_MAX] = {
    0,
    1, 1,
    2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4, 4 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    k = 0;
    x = 0.0;
    b = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    k = k_vec[n_data-1];
    x = x_vec[n_data-1];
    b = b_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *bernstein_poly_ab ( int n, double a, double b, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_AB evaluates at X the Bernstein polynomials based in [A,B].
//
//  Discussion:
//
//    The formula is:
//
//      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
//
//  First values:
//
//    B(0,0)(X) =   1
//
//    B(1,0)(X) = (      B-X                ) / (B-A)
//    B(1,1)(X) = (                 X-A     ) / (B-A)
//
//    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
//    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
//    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
//
//    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
//    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
//    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
//    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
//
//    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
//    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
//    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
//    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
//    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials 
//    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
//    each of degree N, which form a basis for polynomials on [A,B].
//
//    Input, double A, B, the endpoints of the interval on which the
//    polynomials are to be based.  A and B should not be equal.
//
//    Input, double X, the point at which the polynomials 
//    are to be evaluated.
//
//    Output, double BERNSTEIN_POLY_AB[N+1], the values of the N+1
//    Bernstein polynomials at X.
//
{
  double *bern;
  int i;
  int j;

  if ( b == a )
  {
    cerr << "\n";
    cerr << "BERNSTEIN_POLY_AB - Fatal error!\n";
    cerr << "  A = B = " << a << "\n";
    exit ( 1 );
  }

  bern = new double[n+1];

  if ( n == 0 )
  {
    bern[0] = 1.0;
   }
  else if ( 0 < n )
  {
    bern[0] = ( b - x ) / ( b - a );
    bern[1] = ( x - a ) / ( b - a );
 
    for ( i = 2; i <= n; i++ )
    {
      bern[i] = ( x - a ) * bern[i-1] / ( b - a );
      for ( j = i - 1; 1 <= j; j-- )
      {
        bern[j] = ( ( b - x     ) * bern[j]
                  + (     x - a ) * bern[j-1] ) 
                  / ( b     - a );
      }
      bern[0] = ( b - x ) * bern[0] / ( b - a );
    }
  }
  return bern;
}
//****************************************************************************80

double *bernstein_poly_ab_approx ( int n, double a, double b, double ydata[], 
  int nval, double xval[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_AB_APPROX: Bernstein approximant to F(X) on [A,B].
//
//  Formula:
//
//    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
//
//    where
//
//      X(I) = ( ( N - I ) * A + I * B ) / N
//      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
//
//  Discussion:
//
//    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
//    not an interpolant; in other words, its value is not guaranteed to equal
//    that of F at any particular point.  However, for a fixed interval
//    [A,B], if we let N increase, the Bernstein polynomial converges
//    uniformly to F everywhere in [A,B], provided only that F is continuous.
//    Even if F is not continuous, but is bounded, the polynomial converges
//    pointwise to F(X) at all points of continuity.  On the other hand,
//    the convergence is quite slow compared to other interpolation
//    and approximation schemes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomial
//    to be used.  N must be at least 0.
//
//    Input, double A, B, the endpoints of the interval on which the
//    approximant is based.  A and B should not be equal.
//
//    Input, double YDATA[N+1], the data values at N+1 equally
//    spaced points in [A,B].  If N = 0, then the evaluation point should
//    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
//    ( (N-I)*A + I*B ) / N ).
//
//    Input, int NVAL, the number of points at which the
//    approximant is to be evaluated.
//
//    Input, double XVAL[NVAL], the point at which the Bernstein 
//    polynomial approximant is to be evaluated.  The entries of XVAL do not 
//    have to lie in the interval [A,B].
//
//    Output, double BPAB_APPROX[NVAL], the values of the Bernstein 
//    polynomial approximant for F, based in [A,B], evaluated at XVAL.
//
{
  double *bvec;
  int i;
  double *yval;

  yval = new double[nval];

  for ( i = 0; i < nval; i++ )
  {
//
//  Evaluate the Bernstein basis polynomials at XVAL.
//
    bvec = bernstein_poly_ab ( n, a, b, xval[i] );
//
//  Now compute the sum of YDATA(I) * BVEC(I).
//
    yval[i] = r8vec_dot_product ( n + 1, ydata, bvec );

    delete [] bvec;
  }

  return yval;
}
//****************************************************************************80

double *bernstein_vandermonde ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_VANDERMONDE returns the Bernstein Vandermonde matrix.
//
//  Discussion:
//
//    The Bernstein Vandermonde matrix of order N is constructed by
//    evaluating the N Bernstein polynomials of degree N-1 at N equally
//    spaced points between 0 and 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 December 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_VANDERMONDE[N*N], the Bernstein Vandermonde matrix.
//
{
  double *b;
  int i;
  int j;
  double *v;
  double x;

  v = new double[n*n];

  if ( n == 1 )
  {
    v[0+0*1] = 1.0;
    return v;
  }

  for ( i = 0; i < n; i++ )
  {
    x = double( i ) / double( n - 1 );
    b = bernstein_poly_01 ( n - 1, x );
    for ( j = 0; j < n; j++ )
    {
      v[i+j*n] = b[j];
    }
    delete [] b;
  }

  return v;
}
//****************************************************************************80

double r8mat_is_identity ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
//
//  Discussion:
//
//    An R8MAT is a matrix of real ( kind = 8 ) values.
//
//    The routine returns the Frobenius norm of A - I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
//    of the difference matrix A - I, which would be exactly zero
//    if A were the identity matrix.
//
{
  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        t = a[i+j*n] - 1.0;
      }
      else
      {
        t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t * t;
    }
  }
  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
}
//****************************************************************************80

double r8mat_norm_fro ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
