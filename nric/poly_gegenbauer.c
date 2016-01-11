
#include <math.h>

#include "nric.h"


/*******************************************************************************

        gegenbauer_poly

    This routine returns the Gegenbauer or ultraspherical polynomial 
    of order n, alpha: C_n^alpha.

*******************************************************************************/

double  gegenbauer_poly( int n, double alpha, double x)

{
    double  c_0, c_1, c_n, c_n1, c_n2;
    int  nn;

    if  ( n < 0)    nrerror( "Bad argument n in routine gegenbauer_poly()");

    /*
     *    Compute c_0.
     */
    c_0 = 1.0;
    if  ( n == 0)    return  c_0;

    /*
     *    Compute c_1.
     */
    c_1 = 2.0*alpha*x;
    if  ( n == 1)    return  c_1;

    /*
     *    Compute c_n using recurrence on n.
     */
    for  ( c_n2 = c_0, c_n1 = c_1, nn = 2; nn <= n; nn++)    {

        c_n = (1.0/n)*(2.0*(nn - 1.0 + alpha)*x*c_n1 - (1.0*nn - 2.0 + 2.0*alpha)*c_n2);
        c_n2 = c_n1;
        c_n1 = c_n;
    }
    return  c_n;
}

