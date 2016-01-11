
#include <math.h>


#include "nric.h"




double  ln_gamma( double x)  {

    double  xx, temp, series;
    static double  coeff[6] = { 76.18009173,
                               -86.50532033,
                                24.01409822,
                                -1.231739516,
                                 0.120858003e-2,
                                -0.536382e-5};
    int  j;

    xx = x - 1.0;
    temp = xx + 5.5;
    temp -= ( xx + 0.5)*log( temp);
    for  ( series = 1.0, j = 0; j <= 5; j++)    {

        xx += 1.0;
        series += coeff[j]/xx;
    }
    return  -temp + log( 2.50662827465*series);
}



double  factorial( int n)  {

    static int ntop=10;
    static int  a[33] = { 1.0,
                          1.0,
                          2.0,
                          6.0,
                         24.0,
                        120.0,
                        720.0,
                       5040.0,
                      40320.0,
                     362880.0,
                    3628800.0};
    int  j;

    if  ( n < 0)
      nrerror("Negative factorial in routine factorial()");

    if  ( n > 32)
      return  exp( ln_gamma( n + 1.0));

    while  ( ntop < n )  {

        j = ntop++;
        a[ntop] = a[j]*ntop;
    }
    return  a[n];
}



double  bin_coeff( int n, int k)  {

    return  floor( 0.5*exp( ln_factorial( n ) - ln_factorial( k)
                           - ln_factorial( n - k)));
}



double  ln_factorial( int n)  {

    static double  a[101];

    if  ( n < 0 )
      nrerror("Negative factorial in routine ln_factorial().");

    if  ( n <= 1 )
      return  0.0;

    if  ( n <= 100 )
      return  a[n] ? a[n] : ( a[n] = ln_gamma( n + 1.0));

    else
      return  ln_gamma( n + 1.0);
}



double  beta( double z, double w)  {

    return  exp( ln_gamma( z) + ln_gamma( w) - ln_gamma( z + w));
}



double  gamma_p( double a, double x)  {

    double  gamser, gamconfrac, lngam;

    if  ( x < 0.0 || a <= 0.0)
      nrerror( "Invalid arguments in routine gamma_p()");

    if  ( x < ( a + 1.0))    {

        gamma_series( &gamser, a, x, &lngam);
        return  gamser;

    }  else  {

        gamma_cont_fraction( &gamconfrac, a, x, &lngam);
        return  1.0 - gamconfrac;

    }
}



double  gamma_q( double a, double x)

{
    double  gamser, gamconfrac, lngam;

    if  ( x < 0.0 || a <= 0.0)
      nrerror( "Invalid arguments in routine gamma_q().");

    if  ( x < ( a + 1.0))    {

        gamma_series( &gamser, a, x, &lngam);
        return  1.0 - gamser;

    }  else  {

        gamma_cont_fraction( &gamconfrac, a, x, &lngam);
        return  gamconfrac;

    }
}



void  gamma_series( double *gamser, double a, double x, double *lngam)

{
    int  n;
    double  sum, del, ap = 0.0;
    const double EPS = 3.0e-7;
    const int ITMAX = 100;

    *lngam = ln_gamma( a);

    if  ( x <= 0.0)    {

        if  ( x < 0.0)
          nrerror( "Argument less than 0 in routine gamma_series().");
        *gamser = 0.0;
        return;

    }
    else    {

        ap = a;
        del = sum = 1.0/a;
        for  ( n = 1; n <= ITMAX; n++)    {

            ap += 1.0;
            del *= x/ap;
            sum += del;
            if  ( fabs( del) < EPS*fabs( sum))    {

                *gamser = sum*exp( -x + a*log( x) - (*lngam));
                return;
            }
        }
        nrerror( " a too large, ITMAX too small in routine gamma_series().");
        return;
    }
}



void  gamma_cont_fraction( double *gamconfrac,
                           double a, double x, double *lngam)  {

    int  n;
    double  gold = 0.0, g, fact = 1.0, b1 = 1.0;
    double  b0 = 0.0, anf, ana, an, a1 = 0.0, a0 = 1.0;
    const double EPS = 3.0e-7;
    const int ITMAX = 100;

    *lngam = ln_gamma( a);
    a1 = x;
    for  ( n = 1; n <= ITMAX; n++)    {

        an = ( double) n;
        ana = an - a;
        a0 = ( a1 + a0*ana)*fact;
        b0 = ( b1 + b0*ana)*fact;
        anf = an*fact;
        a1 = x*a0 + anf*a1;
        b1 = x*b0 + anf*b1;
        if  ( a1)    {

            fact = 1.0/a1;
            g = b1*fact;
            if  ( fabs( g - gold)/g < EPS)    {

                *gamconfrac = exp( -x + a*log( x) - (*lngam))*g;
                return;
            }
            gold = g;
        }
    }
    nrerror( " a too large, ITMAX too small in routine gamma_series().");
}


double  error_func( double x)  {

    return  ( x < 0) ? ( -gamma_p( 0.5, x*x)) : ( gamma_p( 0.5, x*x));
}


double  error_func_comp( double x)  {

    return  ( x < 0) ? ( 1.0 + gamma_p( 0.5, x*x)) : ( gamma_q( 0.5, x*x));
}
