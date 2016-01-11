
#include <math.h>


#include "nric.h"



/*******************************************************************************

    This routine returns the spherical harmonic of order l, m.

*******************************************************************************/

double  spherical_harmonic( int l, int m, double theta, double phi)

{
    double  Ry_lm, Iy_lm, norm, phase, y_lm;
    int  i, mm, factor1, factor2;


    if  ( ( m > l) || ( m < -l))    nrerror( "Bad arguments in routine legendre_poly()");

    if  ( m < 0)    mm = -m;
    else    mm = m;

/*    Compute modulus of y_lm.    */

    for  ( factor1 = 1, i = 1; i <= ( l - mm); i++)    factor1 *= i;

    for  ( factor2 = 1, i = 1; i <= ( l + mm); i++)    factor2 *= i;

    norm = sqrt( ( ( 2.0*l + 1.0)/ 4.0*PI)*( 1.0*factor1/factor2));

    y_lm = norm*legendre_poly( l, mm, cos( theta));


/*    Compute y_lm for m >= 0.    */

    if  ( m >= 0)    {

        Ry_lm =  y_lm*cos( mm*phi);
        Iy_lm =  y_lm*sin( mm*phi);
    }


/*    Compute y_lm for m < 0.    */

    if  ( m < 0)    {

        for  ( i = 1; i <= mm; i++)    y_lm *= -1.0;
        Ry_lm =  y_lm*cos( mm*phi);
        Iy_lm = -y_lm*sin( mm*phi);
    }

    return  y_lm;
}





/*****************************************************************************************************

    This routine returns the associated Legendre polynomial function of order l, m.
    Here l, m are integers satisfying 0 <= m <= l and x lies in the domain -1.0 <= x <= 1.0.

*****************************************************************************************************/

double  legendre_poly( int l, int m, double x)

{
    double  factor, p_ll, p_mm, p_mmp1, somx2;
    int  i, ll;

    if  ( ( m < 0 ) || ( m > l ) )    nrerror( "Bad argument m in routine legendre_poly()");
    if  ( fabs( x ) > 1.0 )    nrerror( "Bad argument x in routine legendre_poly()");


    /*    Compute P^m_m    */

    p_mm = 1.0;
    if  ( m > 0 )    {
        somx2 = sqrt( ( 1.0 - x )*( 1.0 + x ) );
        for ( factor = 1.0, i = 1; i <= m; i++ )    {
            p_mm *= -factor*somx2;
            factor += 2.0;
        }
    }
    if  ( l == m )    return  p_mm;

    else {

        /*    Compute P^m_m+1    */
        p_mmp1 = ( 2*m + 1 )*x*p_mm;
        if  ( l == ( m + 1 ) )    return  p_mmp1;

        /*    Compute P^m_l, l > m+1    */
        else {
            for  ( ll = ( m + 2 ); ll <= l; ll++ )    {
                p_ll = ( ( 2*ll - 1 )*x*p_mmp1 - ( ll + m - 1 )*p_mm )/( ll - m );
                p_mm = p_mmp1;
                p_mmp1 = p_ll;
            }
            return  p_ll;
        }
    }
}

