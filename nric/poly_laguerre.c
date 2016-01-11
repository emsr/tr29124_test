
#include <math.h>

#include "nric.h"


/*******************************************************************************

        LAGUERRE_POLY

    This routine returns the associated Laguerre polynomial 
    of order n, alpha: L_n^alpha.

*******************************************************************************/

double  laguerre_poly( int n, double alpha, double x)

{
    double  l_0, l_1, l_n, l_n1, l_n2;
    int  nn;

    if  ( n < 0)    nrerror( "Bad argument n in routine laguerre_poly()");

    /*
     *    Compute l_0.
     */
    l_0 = 1.0;
    if  ( n == 0)    return  l_0;


    /*
     *    Compute l_1.
     */
    l_1 = -x + 1.0 + alpha;
    if  ( n == 1)    return  l_1;


    /*
     *    Compute l_n by recursion on n.
     */
    for  ( l_n2 = l_0, l_n1 = l_1, nn = 2; nn <= n; nn++)    {

        l_n = ( 2.0*nn - 1.0 + alpha - x)*l_n1 - ( 1.0*nn - 1.0 + alpha)*l_n2;
        l_n2 = l_n1;
        l_n1 = l_n;
    }
    return  l_n;
}

