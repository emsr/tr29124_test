# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "gegenbauer_polynomial.hpp"

//****************************************************************************80

bool gegenbauer_alpha_check ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_ALPHA_CHECK checks the value of ALPHA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, a parameter which is part of the definition of
//    the Gegenbauer polynomials.  It must be greater than -0.5.
//
//    Output, bool GEGENBAUER_ALPHA_CHECK.
//    TRUE, ALPHA is acceptable.
//    FALSE, ALPHA is not acceptable. 
//
{
  bool check;
  bool squawk;

  squawk = false;

  if ( -0.5 < alpha )
  {
    check = true;
  }
  else
  {
    check = false;
    if ( squawk )
    {
      cerr << "\n";
      cerr << "GEGENBAUER_ALPHA_CHECK - Fatal error!\n";
      cerr << "  Illegal value of ALPHA.\n";
      cerr << "  ALPHA = " << alpha << "\n";
      cerr << "  but ALPHA must be greater than -0.5.\n";
    }
  }

  return check;
}
//****************************************************************************80

void gegenbauer_ek_compute ( int n, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_EK_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2015
//
//  Author:
//
//    John Burkardt
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
//    Input, int N, the order of the quadrature rule.
//
//    Input, double ALPHA, the exponent of (1-X^2) in the weight.  
//    -1.0 < ALPHA is required.
//
//    Input, double A, B, the left and right endpoints 
//    of the interval.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double abi;
  double *bj;
  bool check;
  int i;
  double zemu;
//
//  Check N.
//
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_EK_COMPUTE - Fatal error!\n";
    cerr << "  1 <= N is required.\n";
    exit ( 1 );
  }
//
//  Check ALPHA.
//
  check = gegenbauer_alpha_check ( alpha );
  if ( ! check )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_EK_COMPUTE - Fatal error!\n";
    cerr << "  Illegal value of ALPHA.\n";
    exit ( 1 );
  }
//
//  Define the zero-th moment.
//
  zemu = pow ( 2.0, 2.0 * alpha + 1.0 )
    * r8_gamma ( alpha + 1.0 )
    * r8_gamma ( alpha + 1.0 )
    / r8_gamma ( 2.0 * alpha + 2.0 );
//
//  Define the Jacobi matrix.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  bj = new double[n];

  bj[0] = 4.0 * pow ( alpha + 1.0, 2 )
    / ( ( 2.0 * alpha + 3.0 ) * pow ( 2.0 * alpha + 2.0, 2 ) );

  for ( i = 2; i <= n; i++ )
  {
    abi = 2.0 * ( alpha + ( double ) i );
    bj[i-1] = 4.0 * ( double ) ( i ) * pow ( alpha + i, 2 ) * ( 2.0 * alpha + i )
      / ( ( abi - 1.0 ) * ( abi + 1.0 ) * abi * abi );
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
    w[i] = pow ( w[i], 2 );
  }

  free ( bj );

  return;
}
//****************************************************************************80

double gegenbauer_integral ( int expon, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_INTEGRAL evaluates the integral of a monomial with Gegenbauer weight.
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x <= +1 ) x^expon (1-x^2)^alpha dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON, the exponent.
//
//    Input, double ALPHA, the exponent of (1-X^2) in the weight factor.
//
//    Output, double GEGENBAUER_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = tgamma ( 1.0 + c ) * 2.0 
    * tgamma ( 1.0 + alpha  ) * value1 
    / tgamma ( 2.0 + alpha  + c );

  return value;
}
//****************************************************************************80

double *gegenbauer_polynomial_value ( int m, int n, double alpha, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLYNOMIAL_VALUE computes the Gegenbauer polynomials C(I,ALPHA)(X).
//
//  Differential equation:
//
//    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + M (M + 2 ALPHA) Y = 0
//
//  Recursion:
//
//    C(0,ALPHA,X) = 1,
//    C(1,ALPHA,X) = 2*ALPHA*X
//    C(M,ALPHA,X) = (  ( 2*M-2+2*ALPHA) * X * C(M-1,ALPHA,X) 
//                    + (  -M+2-2*ALPHA)   *   C(M-2,ALPHA,X) ) / M
//
//  Restrictions:
//
//    ALPHA must be greater than -0.5.
//
//  Special values:
//
//    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
//    polynomials of the second kind.
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X^2 )^( ALPHA - 0.5 ) * C(M,ALPHA,X)^2 dX
//
//    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( M + 2 * ALPHA ) 
//      / ( M! * ( M + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2015
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
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input, int M, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double ALPHA, a parameter which is part of the definition of
//    the Gegenbauer polynomials.  It must be greater than -0.5.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double GEGENBAUER_POLYNOMIAL_VALUE(1:M+1,N), the values of 
//    Gegenbauer polynomials 0 through M
//    at the N points X.  
//
{
  double *c;
  bool check;
  int i;
  double i_r8;
  int j;

  check = gegenbauer_alpha_check ( alpha );
  if ( ! check )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_POLYNOMIAL_VALUE - Fatal error!\n";
    cerr << "  Illegal value of ALPHA.\n";
    exit ( 1 );
  }

  c = new double[(m+1)*n];

  if ( m < 0 )
  {
    return c;
  }

  if ( n == 0 )
  {
    return c;
  }

  for ( j = 0; j < n; j++ )
  {
    c[0+j*(m+1)] = 1.0;
  }

  if ( m < 1 )
  {
    return c;
  }

  for ( j = 0; j < n; j++ )
  {
    c[1+j*(m+1)] = 2.0 * alpha * x[j];
  }

  for ( i = 2; i <= m; i++ )
  {
    i_r8 = ( double ) i;
    for ( j = 0; j < n; j++ )
    {
      c[i+j*(m+1)] = (  (     2.0 * i_r8 - 2.0  + 2.0 * alpha ) * x[j] * c[i-1+j*(m+1)]
                     +  (         - i_r8 + 2.0  - 2.0 * alpha ) *        c[i-2+j*(m+1)] )
                     /              i_r8 ;
    }
  }

  return c;
}
//****************************************************************************80

void gegenbauer_polynomial_values ( int &n_data, int &n, double &a, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLYNOMIAL_VALUES returns some values of the Gegenbauer polynomials.
//
//  Discussion:
//
//    The Gegenbauer polynomials are also known as the "spherical
//    polynomials" or "ultraspherical polynomials".
//
//    In Mathematica, the function can be evaluated by:
//
//      GegenbauerC[n,m,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2004
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
//    Output, int &N, the order parameter of the function.
//
//    Output, double &A, the real parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 38

  static double a_vec[N_MAX] = {
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
      0.0E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00 };

  static double fx_vec[N_MAX] = {
       1.0000000000E+00,
       0.2000000000E+00,
      -0.4400000000E+00,
      -0.2800000000E+00,
       0.2320000000E+00,
       0.3075200000E+00,
      -0.0805760000E+00,
      -0.2935168000E+00,
      -0.0395648000E+00,
       0.2459712000E+00,
       0.1290720256E+00,
       0.0000000000E+00,
      -0.3600000000E+00,
      -0.0800000000E+00,
       0.8400000000E+00,
       2.4000000000E+00,
       4.6000000000E+00,
       7.4400000000E+00,
      10.9200000000E+00,
      15.0400000000E+00,
      19.8000000000E+00,
      25.2000000000E+00,
      -9.0000000000E+00,
      -0.1612800000E+00,
      -6.6729600000E+00,
      -8.3750400000E+00,
      -5.5267200000E+00,
       0.0000000000E+00,
       5.5267200000E+00,
       8.3750400000E+00,
       6.6729600000E+00,
       0.1612800000E+00,
      -9.0000000000E+00,
     -15.4252800000E+00,
      -9.6969600000E+00,
      22.4409600000E+00,
     100.8892800000E+00,
     252.0000000000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  2,
     2,  2,  2,
     2,  2,  2,
     2,  2,  2,
     2,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5 };

  static double x_vec[N_MAX] = {
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
     -0.50E+00,
     -0.40E+00,
     -0.30E+00,
     -0.20E+00,
     -0.10E+00,
      0.00E+00,
      0.10E+00,
      0.20E+00,
      0.30E+00,
      0.40E+00,
      0.50E+00,
      0.60E+00,
      0.70E+00,
      0.80E+00,
      0.90E+00,
      1.00E+00 };

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
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gegenbauer_ss_compute ( int order, double alpha, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_SS_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//    Thanks to Janiki Raman for pointing out a problem in an earlier
//    version of the code that occurred when ALPHA was -0.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule.
//
//    Input, double ALPHA, the exponent of (1-X^2) in the weight.  
//    -1.0 < ALPHA is required.
//
//    Output, double XTAB[ORDER], the abscissas.
//
//    Output, double WEIGHT[ORDER], the weights.
//
{
  double an;
  double *c;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x;
//
//  Check ORDER.
//
  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_SS_COMPUTE - Fatal error!\n";
    cerr << "  1 <= ORDER is required.\n";
    exit ( 1 );
  }
  
  c = new double[order];
//
//  Check ALPHA.
//
  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_SS_COMPUTE - Fatal error!\n";
    cerr << "  -1.0 < ALPHA is required.\n";
    exit ( 1 );
  }
//
//  Set the recursion coefficients.
//
  c[0] = 0.0;
  if ( 2 <= order )
  {
    c[1] = 1.0 / ( 2.0 * alpha + 3.0 );
  }

  for ( i = 3; i <= order; i++ )
  {
    c[i-1] = ( double ) ( i - 1 ) 
          * ( alpha + alpha + ( double ) ( i - 1 ) ) / 
          ( ( alpha + alpha + ( double ) ( 2 * i - 1 ) ) 
          * ( alpha + alpha + ( double ) ( 2 * i - 3 ) ) );
  }

  delta = tgamma ( alpha         + 1.0 ) 
        * tgamma (         alpha + 1.0 ) 
        / tgamma ( alpha + alpha + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + alpha + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );

      r1 = ( 1.0 + alpha ) 
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) ) 
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 2.44 * an + 1.282 * an * an;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) / 
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) * 
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * alpha * 
        ( 1.0 + 0.25 * fabs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * alpha ) / ( 0.766 + 0.119 * alpha );

      r2 = 1.0 / ( 1.0 + 0.639 
        * ( ( double ) ( order ) - 4.0 ) 
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) * 
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * alpha ) / ( 1.67 + 0.28 * alpha );

      r2 = 1.0 / 
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 ) 
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha / 
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }

    gegenbauer_ss_root ( x, order, alpha, dp2, p1, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
//
//  Reverse the order of the values.
//
  for ( i = 1; i <= order/2; i++ )
  {
    temp          = xtab[i-1];
    xtab[i-1]     = xtab[order-i];
    xtab[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp            = weight[i-1];
    weight[i-1]     = weight[order-i];
    weight[order-i] = temp;
  }

  delete [] c;

  return;
}
//****************************************************************************80

void gegenbauer_ss_recur ( double &p2, double &dp2, double &p1, double x,
  int order, double alpha, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_SS_RECUR: value and derivative of a Gegenbauer polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double &P2, the value of J(ORDER)(X).
//
//    Output, double &DP2, the value of J'(ORDER)(X).
//
//    Output, double &P1, the value of J(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponents of (1-X^2).
//
//    Input, double C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  p1 = 1.0;
  dp1 = 0.0;

  p2 = x;
  dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = p1;
    dp0 = dp1;

    p1 = p2;
    dp1 = dp2;

    p2 = x *  p1 - c[i-1] * p0;
    dp2 = x * dp1 + p1 - c[i-1] * dp0;
  }
  return;
}
//****************************************************************************80

void gegenbauer_ss_root ( double &x, int order, double alpha,  double &dp2, 
  double &p1, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_SS_ROOT improves a root of a Gegenbauer polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double &X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponents of (1-X^2).
//
//    Output, double &DP2, the value of J'(ORDER)(X).
//
//    Output, double &P1, the value of J(ORDER-1)(X).
//
//    Input, double C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gegenbauer_ss_recur ( p2, dp2, p1, x, order, alpha, c );

    d = p2 / dp2;
    x = x - d;

    if ( fabs ( d ) <= eps * ( fabs ( x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
//****************************************************************************80

void hyper_2f1_values ( int &n_data, double &a, double &b, double &c,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = Hypergeometric2F1 [ a, b, c, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2007
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
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, &C, &X, the parameters of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 24

  static double a_vec[N_MAX] = {
   -2.5,
   -0.5,
    0.5,
    2.5,
   -2.5,
   -0.5,
    0.5,
    2.5,
   -2.5,
   -0.5,
    0.5,
    2.5,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3 };
  static double b_vec[N_MAX] = {
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7 };
  static double c_vec[N_MAX] = {
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
   -5.5,
   -0.5,
    0.5,
    4.5,
   -5.5,
   -0.5,
    0.5,
    4.5,
   -5.5,
   -0.5,
    0.5,
    4.5 };
  static double fx_vec[N_MAX] = {
    0.72356129348997784913,
    0.97911109345277961340,
    1.0216578140088564160,
    1.4051563200112126405,
    0.46961431639821611095,
    0.95296194977446325454,
    1.0512814213947987916,
    2.3999062904777858999,
    0.29106095928414718320,
    0.92536967910373175753,
    1.0865504094806997287,
    5.7381565526189046578,
    15090.669748704606754,
   -104.31170067364349677,
    21.175050707768812938,
    4.1946915819031922850,
    1.0170777974048815592E+10,
   -24708.635322489155868,
    1372.2304548384989560,
    58.092728706394652211,
    5.8682087615124176162E+18,
   -4.4635010147295996680E+08,
    5.3835057561295731310E+06,
    20396.913776019659426 };
  static double x_vec[N_MAX] = {
    0.25,
    0.25,
    0.25,
    0.25,
    0.55,
    0.55,
    0.55,
    0.55,
    0.85,
    0.85,
    0.85,
    0.85,
    0.25,
    0.25,
    0.25,
    0.25,
    0.55,
    0.55,
    0.55,
    0.55,
    0.85,
    0.85,
    0.85,
    0.85 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    c = c_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
