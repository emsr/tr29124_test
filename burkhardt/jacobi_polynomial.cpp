# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <sstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "jacobi_polynomial.hpp"

//****************************************************************************80

string i4_to_string ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}
//****************************************************************************80

double j_double_product_integral ( int i, int j, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    J_DOUBLE_PRODUCT_INTEGRAL: integral of J(i,x)*J(j,x)*(1-x)^a*(1+x)^b.
//
//  Discussion:
//
//    VALUE = integral ( -1 <= x <= +1 ) J(i,x)*J(j,x)*(1-x)^a*(1+x)^b dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//
//    Input, double A, B, the parameters.
//    -1 < A, B.
//
//    Output, double VALUE, the value of the integral.
//
{ 
  double i_r8;
  double value;

  if ( i != j )
  {
    value = 0.0;
  }
  else
  {
    i_r8 = ( double ) ( i );

    value = pow ( 2, a + b + 1.0 ) 
      / ( 2.0 * i_r8 + a + b + 1.0 ) 
      * tgamma ( i_r8 + a + 1.0 ) 
      * tgamma ( i_r8 + b + 1.0 )
      / r8_factorial ( i ) 
      / tgamma ( i_r8 + a + b + 1.0 );
  }
  return value;
}
//****************************************************************************80

double j_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    J_INTEGRAL evaluates a monomial integral associated with J(n,a,b,x).
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x < +1 ) x^n dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the exponent.
//    0 <= N.
//
//    Output, double J_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 / ( double ) ( n + 1 );
  }

  return value;
}
//****************************************************************************80

double *j_polynomial ( int m, int n, double alpha, double beta, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY evaluates the Jacobi polynomial J(n,a,b,x).
//
//  Differential equation:
//
//    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
//
//  Recursion:
//
//    P(0,ALPHA,BETA,X) = 1,
//
//    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
//
//    P(N,ALPHA,BETA,X)  =
//      (
//        (2*N+ALPHA+BETA-1)
//        * ((ALPHA^2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
//        * P(N-1,ALPHA,BETA,X)
//        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
//      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
//
//  Restrictions:
//
//    -1 < ALPHA
//    -1 < BETA
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
//      * P(N,ALPHA,BETA,X)^2 dX
//    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
//      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
//
//  Special values:
//
//    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.  Note
//    that polynomials 0 through N will be computed.
//
//    Input, double ALPHA, one of the parameters defining the Jacobi
//    polynomials, ALPHA must be greater than -1.
//
//    Input, double BETA, the second parameter defining the Jacobi
//    polynomials, BETA must be greater than -1.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double J_POLYNOMIAL[M*(N+1)], the values.
//
{
  double c1;
  double c2;
  double c3;
  double c4;
  int i;
  int j;
  double *v;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "J_POLYNOMIAL - Fatal error!\n";
    cerr << "  Illegal input value of ALPHA = " << alpha << "\n";
    cerr << "  But ALPHA must be greater than -1.\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "J_POLYNOMIAL - Fatal error!\n";
    cerr << "  Illegal input value of BETA = " << beta << "\n";
    cerr << "  But BETA must be greater than -1.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = ( 1.0 + 0.5 * ( alpha + beta ) ) * x[i] 
      + 0.5 * ( alpha - beta );
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      c1 = 2.0 * ( double ) ( j ) * ( ( double ) ( j ) + alpha + beta ) 
        * ( ( double ) ( 2 * j - 2 ) + alpha + beta );

      c2 = ( ( double ) ( 2 * j - 1 ) + alpha + beta ) 
        * ( ( double ) ( 2 * j ) + alpha + beta ) 
        * ( ( double ) ( 2 * j - 2 ) + alpha + beta );

      c3 = ( ( double ) ( 2 * j - 1 ) + alpha + beta ) 
        * ( alpha + beta ) * ( alpha - beta );

      c4 = - ( double ) ( 2 ) * ( ( double ) ( j - 1 ) + alpha ) 
        * ( ( double ) ( j - 1 ) + beta )  
        * ( ( double ) ( 2 * j ) + alpha + beta );

      v[i+j*m] = ( ( c3 + c2 * x[i] ) * v[i+(j-1)*m] + c4 * v[i+(j-2)*m] ) / c1;
    }
  }

  return v;
}
//****************************************************************************80

void j_polynomial_values ( int &n_data, int &n, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    J_POLYNOMIAL_VALUES returns some values of the Jacobi polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiP[ n, a, b, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
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
//    Output, double &A, &B, parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 26

  static double a_vec[N_MAX] = {
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 1.0, 2,
     3.0, 4.0, 5.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0, 0.0, 0,
     0.0, 0.0 };

  static double b_vec[N_MAX] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 2.0,
    3.0, 4.0, 5.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.3750000000000000E+00,
     -0.4843750000000000E+00,
     -0.1328125000000000E+00,
      0.2753906250000000E+00,
     -0.1640625000000000E+00,
     -0.1174804687500000E+01,
     -0.2361328125000000E+01,
     -0.2616210937500000E+01,
      0.1171875000000000E+00,
      0.4218750000000000E+00,
      0.5048828125000000E+00,
      0.5097656250000000E+00,
      0.4306640625000000E+00,
     -0.6000000000000000E+01,
      0.3862000000000000E-01,
      0.8118400000000000E+00,
      0.3666000000000000E-01,
     -0.4851200000000000E+00,
     -0.3125000000000000E+00,
      0.1891200000000000E+00,
      0.4023400000000000E+00,
      0.1216000000000000E-01,
     -0.4396200000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0, 1, 2, 3,
     4, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
     -1.0E+00,
     -0.8E+00,
     -0.6E+00,
     -0.4E+00,
     -0.2E+00,
      0.0E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *j_polynomial_zeros ( int n, double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int, N, the order.
//
//    Input, double, ALPHA, BETA, the parameters.
//    -1 < ALPHA, BETA.
//
//    Output, double J_POLYNOMIAL_ZEROS[N], the zeros.
//
{
  double a2b2;
  double ab;
  double abi;
  double *bj;
  int i;
  double i_r8;
  double *w;
  double *x;
  double zemu;

  ab = alpha + beta;
  abi = 2.0 + ab;
//
//  Define the zero-th moment.
//
  zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 ) 
    * tgamma ( beta + 1.0 ) / tgamma ( abi );
//
//  Define the Jacobi matrix.
//
  x = new double[n];
  x[0] = ( beta - alpha ) / abi;
  for ( i = 1; i < n; i++ )
  {
    x[i] = 0.0;
  }

  bj = new double[n];

  bj[0] = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
    / ( ( abi + 1.0 ) * abi * abi );
  for ( i = 1; i < n; i++ )
  {
    bj[i] = 0.0;
  }

  a2b2 = beta * beta - alpha * alpha;

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    abi = 2.0 * i_r8 + ab;
    x[i] = a2b2 / ( ( abi - 2.0 ) * abi );
    abi = abi * abi;
    bj[i] = 4.0 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) 
      * ( i_r8 + ab ) / ( ( abi - 1.0 ) * abi );
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
  }

  w = new double[n];

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  delete [] bj;
  delete [] w;

  return x;
}
//****************************************************************************80

void j_quadrature_rule ( int n, double alpha, double beta, double x[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int, N, the order.
//
//    Input, double, ALPHA, BETA, the parameters.
//    -1 < ALPHA, BETA.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double a2b2;
  double ab;
  double abi;
  double *bj;
  int i;
  double i_r8;
  double zemu;

  ab = alpha + beta;
  abi = 2.0 + ab;
//
//  Define the zero-th moment.
//
  zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 ) 
    * tgamma ( beta + 1.0 ) / tgamma ( abi );
//
//  Define the Jacobi matrix.
//
  x[0] = ( beta - alpha ) / abi;
  for ( i = 1; i < n; i++ )
  {
    x[i] = 0.0;
  }

  bj = new double[n];

  bj[0] = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
    / ( ( abi + 1.0 ) * abi * abi );
  for ( i = 1; i < n; i++ )
  {
    bj[i] = 0.0;
  }

  a2b2 = beta * beta - alpha * alpha;

  for ( i = 1; i < n; i++ )
  {
    i_r8 = ( double ) ( i + 1 );
    abi = 2.0 * i_r8 + ab;
    x[i] = a2b2 / ( ( abi - 2.0 ) * abi );
    abi = abi * abi;
    bj[i] = 4.0 * i_r8 * ( i_r8 + alpha ) * ( i_r8 + beta ) 
      * ( i_r8 + ab ) / ( ( abi - 1.0 ) * abi );
  }

  for ( i = 0; i < n; i++ )
  {
    bj[i] = sqrt ( bj[i] );
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
