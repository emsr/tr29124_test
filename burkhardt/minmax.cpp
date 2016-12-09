#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int i1 and i2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
//
{
  if ( i1 > i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int i1 and i2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
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
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  double value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    Factorial2(N)
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial
//    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
//
//    Output, double R8_FACTORIAL2, the value of Factorial2(N).
//
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  const double r8_pi = 3.1415926535897932384626434;
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - r8_pi / sin ( r8_pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8_psi ( double xx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the function Psi(X).
//
//  Discussion:
//
//    This routine evaluates the logarithmic derivative of the
//    Gamma function,
//
//      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
//             = d/dX LN ( GAMMA(X) )
//
//    for real X, where either
//
//      - XMAX1 < X < - XMIN, and X is not a negative integer,
//
//    or
//
//      XMIN < X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody, Anthony Strecok, Henry Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, Number 121, January 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double XX, the argument of the function.
//
//    Output, double R8_PSI, the value of the function.
//
{
  double aug;
  double den;
  double four = 4.0;
  double fourth = 0.25;
  double half = 0.5;
  int i;
  int n;
  int nq;
  double one = 1.0;
  double p1[9] = { 
   4.5104681245762934160E-03, 
   5.4932855833000385356, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = { 
  -2.7103228277757834192, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821, 
  -1.4479614616899842986, 
  -7.3689600332394549911E-02, 
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = { 
   9.6141654774222358525E+01, 
   2.6287715790581193330E+03, 
   2.9862497022250277920E+04, 
   1.6206566091533671639E+05, 
   4.3487880712768329037E+05, 
   5.4256384537269993733E+05, 
   2.4242185002017985252E+05, 
   6.4155223783576225996E-08 };
  double q2[6] = { 
   4.4992760373789365846E+01, 
   2.0240955312679931159E+02, 
   2.4736979003315290057E+02, 
   1.0742543875702278326E+02, 
   1.7463965060678569906E+01, 
   8.8427520398873480342E-01 };
  double sgn;
  double three = 3.0;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;
  double zero = 0.0;

  x = xx;
  w = fabs ( x );
  aug = zero;
//
//  Check for valid arguments, then branch to appropriate algorithm.
//
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( zero < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < half )
  {
//
//  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
//  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
//
    if ( w <= xsmall )
    {
      aug = - one / x;
    }
//
//  Argument reduction for cotangent.
//
    else
    {
      if ( x < zero )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = int ( w * four );
      w = four * ( w - ( double ) ( nq ) * fourth );
//
//  W is now related to the fractional part of 4.0 * X.
//  Adjust argument to correspond to values in the first
//  quadrant and determine the sign.
//
      n = nq / 2;

      if ( n + n != nq )
      {
        w = one - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
//
//  Determine the final value for  -pi * cotan(pi*x).
//
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
//
//  Check for singularity.
//
        if ( z == zero )
        {
          if ( zero < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
          }
          return value;
        }
        aug = sgn * ( four / tan ( z ) );
      }
      else
      {
        aug = sgn * ( four * tan ( z ) );
      }
    }
    x = one - x;
  }
//
//  0.5 <= X <= 3.0.
//
  if ( x <= three )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
//
//  3.0 < X.
//
  if ( x < xlarge )
  {
    w = one / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - half / x + aug;
  }

  value = aug + log ( x );

  return value;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  std::cout << "\n";
  std::cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    std::cout << "\n";
    std::cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    std::cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    std::cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      std::cout << std::setw(7) << j - 1 << "       ";
    }
    std::cout << "\n";
    std::cout << "  Row\n";
    std::cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      std::cout << std::setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        std::cout << std::setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      std::cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
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
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

std::string r8_to_string ( double r8, std::string format )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_STRING converts an R8 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R8, a double.
//
//    Input, string FORMAT, the format string.
//
//    Output, string R8_TO_STRING, the string.
//
{
  char r8_char[80];
  std::string r8_string;

  sprintf ( r8_char, format.c_str ( ), r8 );

  r8_string = std::string ( r8_char );

  return r8_string;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
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
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  std::cout << "\n";
  std::cout << title << "\n";
  std::cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    std::cout << "  " << std::setw(8)  << i
         << ": " << std::setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec2_print ( int n, double a1[], double a2[], std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  std::cout << "\n";
  std::cout << title << "\n";
  std::cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    std::cout << std::setw(6)  << i
         << ": " << std::setw(14) << a1[i]
         << "  " << std::setw(14) << a2[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r8poly_print ( int n, double a[], std::string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT prints out a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X^N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  double mag;
  int n2;
  char plus_minus;

  if ( 0 < title.length ( ) )
  {
    std::cout << "\n";
    std::cout << title << "\n";
  }
  std::cout << "\n";

  if ( n < 0 )
  {
    std::cout << "  p(x) = 0\n";
    return;
  }

  if ( a[n] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = fabs ( a[n] );

  if ( 2 <= n )
  {
    std::cout << "  p(x) = " << plus_minus
         << std::setw(14) << mag << " * x ^ " << n << "\n";
  }
  else if ( n == 1 )
  {
    std::cout << "  p(x) = " << plus_minus
         << std::setw(14) << mag << " * x\n";
  }
  else if ( n == 0 )
  {
    std::cout << "  p(x) = " << plus_minus
         << std::setw(14) << mag << "\n";
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = fabs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        std::cout << "         " << plus_minus
             << std::setw(14) << mag << " * x ^ " << i << "\n";
      }
      else if ( i == 1 )
      {
        std::cout << "         " << plus_minus
             << std::setw(14) << mag << " * x\n";
      }
      else if ( i == 0 )
      {
        std::cout << "         " << plus_minus
             << std::setw(14) << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8mat_mtv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8MAT_MTV_NEW[N], the product A'*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

double r8_hyper_2f1 ( double a, double b, double c, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
//
//  Discussion:
//
//    A minor bug was corrected.  The HW variable, used in several places as
//    the "old" value of a quantity being iteratively improved, was not
//    being initialized.  JVB, 11 February 2008.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    Original FORTRAN77 original by Shanjie Zhang, Jianming Jin.
//    C++ version by John Burkardt.
//
//    The original FORTRAN77 version of this routine is copyrighted by
//    Shanjie Zhang and Jianming Jin.  However, they give permission to
//    incorporate this routine into a user program provided that the copyright
//    is acknowledged.
//
//  Reference:
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//  Parameters:
//
//    Input, double A, B, C, X, the arguments of the function.
//    C must not be equal to a nonpositive integer.
//    X < 1.
//
//    Output, double R8_HYPER_2F1, the value of the function.
//
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  bool l0;
  bool l1;
  bool l2;
  bool l3;
  bool l4;
  bool l5;
  int m;
  int nm;
  double pa;
  double pb;
  const double r8_pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 )
  {
    std::cerr << "\n";
    std::cerr << "R8_HYPER_2F1 - Fatal error!\n";
    std::cerr << "  The hypergeometric series is divergent.\n";
    std::cerr << "  C is integral and negative.\n";
    std::cerr << "  C = " << c << "\n";
    exit ( 1 );
  }

  if ( l1 )
  {
    std::cerr << "\n";
    std::cerr << "R8_HYPER_2F1 - Fatal error!\n";
    std::cerr << "  The hypergeometric series is divergent.\n";
    std::cerr << "  1 - X < 0, C - A - B <= 0\n";
    std::cerr << "  A = " << a << "\n";
    std::cerr << "  B = " << b << "\n";
    std::cerr << "  C = " << c << "\n";
    std::cerr << "  X = " << x << "\n";
    exit ( 1 );
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = tgamma ( c );
    gcab = tgamma ( c - a - b );
    gca = tgamma ( c - a );
    gcb = tgamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && fabs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( r8_pi ) * pow ( 2.0, - a );
    g1 = tgamma ( c );
    g2 = tgamma ( 1.0 + a / 2.0 - b );
    g3 = tgamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( fabs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( fabs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( fabs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( fabs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( fabs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = ( int ) ( c - a - b );
      ga = tgamma ( a );
      gb = tgamma ( b );
      gc = tgamma ( c );
      gam = tgamma ( a + m );
      gbm = tgamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 ) 
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a ) 
              / ( ( j + k ) * ( a + j + k - 1.0 ) ) 
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( fabs ( f1 - hw ) < fabs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) 
            / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( fabs ( f1 - hw ) < fabs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = tgamma ( a );
      gb = tgamma ( b );
      gc = tgamma ( c );
      gca = tgamma ( c - a );
      gcb = tgamma ( c - b );
      gcab = tgamma ( c - a - b );
      gabc = tgamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( fabs ( hf - hw ) < fabs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( fabs ( hf - hw ) <= fabs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    std::cout << "\n";
    std::cout << "R8_HYPER_2F1 - Warning!\n";
    std::cout << "  A large number of iterations were needed.\n";
    std::cout << "  The accuracy of the results should be checked.\n";
  }

  return hf;
}
//****************************************************************************80

double r8_uniform_ab ( double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
  const int i4_huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    std::cerr << "\n";
    std::cerr << "R8_UNIFORM_AB - Fatal error!\n";
    std::cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and the value will be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( ( double ( *seed ) ) * 4.656612875E-10 );

  return r;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    std::cerr << "\n";
    std::cerr << "I4_UNIFORM_AB - Fatal error!\n";
    std::cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentially) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        std::cerr << "\n";
        std::cerr << "IMTQLX - Fatal error!\n";
        std::cerr << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r = sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r = sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r = sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
