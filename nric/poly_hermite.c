
#include <stdio.h>
#include <math.h>

#include "nric.h"

/*******************************************************************************

        HERMITE_POLY

    This routine returns the Hermite polynomial of order n: H_n.

*******************************************************************************/


double  hermite_h0( double x ) {

    return 1.0;
}


double  hermite_h1( double x ) {

    return 2.0*x;
}


double  hermite_h( int n, double x)

{
    int i;
    double  h_0, h_1, h_n, h_n1, h_n2;


    if  ( n < 0)    nrerror( "Bad argument n in routine hermite_poly()");

    /*    Compute h_0    */
    h_0 = 1.0;
    if  ( n == 0)    return h_0;

    /*    Compute h_1    */
    h_1 = 2.0*x;
    if  ( n == 1)    return h_1;

    /*    Compute h_n    */
    for  ( h_n2 = h_0, h_n1 = h_1, i = 2; i <= n; i++)    {

        h_n = 2.0*( x*h_n1 + (i - 1)*h_n2 );
        h_n2 = h_n1;
        h_n1 = h_n;
    }
    return  h_n;
}

