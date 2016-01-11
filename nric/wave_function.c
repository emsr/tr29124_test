
#include <math.h>


#include "nric.h"


/*******************************************************************************

        RADIAL_HYDROGEN

    This routine returns the radial wave function of a one-electron atom: R_nl.
  As a bonus it computes the energy eigenvalue in eV and puts it in a 
  pointer *E_n.

*******************************************************************************/

double  radial_Hydrogen( int Z, int n, int l, double r, double *E_n)

{
    int  factor1, factor2, i;
    double  mu, R_nl, rho, r_Bohr, r_0, norm, phase;
    const double  alpha = 0.30282208, hbarc = 197.329, me = 511000.0;

    if( n <= 0)  nrerror( "Bad radial quantum number n in radial_Hydrogen()");
    if( l >= n)  nrerror( "Bad orbital quantum number l in radial_Hydrogen()");
    if( r < 0.0)  nrerror( "Bad argument r in radial_Hydrogen()");


/*    Calculate Bohr radius.    */

    r_Bohr = hbarc/( me*alpha);
    r_0 = ( n*r_Bohr)/( 2.0*Z);


/*    Calculate energy.    */

    *E_n = -Z*Z*alpha*hbarc/( 2.0*r_Bohr*n*n);


/*    Calculate normalization.    */

    for  ( factor1 = 1, i = 1; i <= ( n - l - 1); i++)    factor1 *= i;
    for  ( factor2 = 1, i = 1; i <= ( n + l); i++)    factor2 *= i;
    norm = sqrt( 1.0*factor1/( 2.0*n*factor2*factor2*factor2*r_0*r_0*r_0));


/*    Define phase    */

    phase = -1.0;


/*    Calculate eigenfunction.    */

    rho = r/r_0;

    R_nl = phase*norm;
    for  ( i = 1; i <= l; i++)    R_nl *= rho;
    R_nl *= exp( rho/2.0);
    R_nl *= laguerre_poly( n + l, ( float) ( 2*l + 1), rho);

    return  R_nl;
}






/*******************************************************************************

        HARMONIC

    This routine returns the harmonic oscillator wave function 
    of order order n: F_n.

*******************************************************************************/

double  harmonic( double mu, double omega, int n, double x, double *E_n)

{
    int  i, factor2;
    double  alpha, norm, F_n;
    const double  hbarc = 197.329;


    if  ( n < 0)    nrerror( "Bad argument n in harmonic().");

    alpha = sqrt( mu*omega)/hbarc;
    if  ( alpha == 0.0 )    return  0.0;


/*    Calculate energy.    */

    *E_n = omega*( n + 0.5);


/*    Calculate normalization.    */

    for  ( factor2 = 1, i = 1; i <= n; i++)    factor2 *= 2*i;
    norm = sqrt( alpha/( factor2*sqrt( PI)));


/*    Calculate wavefunction.    */

    x *= alpha;

    F_n = norm;
    F_n *= exp( -x*x/2.0);
    F_n *= hermite_h( n, x);

    return  F_n;
}






/*******************************************************************************

        RADIAL_HARMONIC_2D

    This routine returns the 2-dimensional harmonic oscillator 
    wave function of order order n: H_n.

*******************************************************************************/

double  radial_harmonic_2d( double mu, double omega, int n, int m, 
                            double p, double *E_nm)

{
    int  i, mm, factor1, factor2;
    double  alpha, x, norm, P_nm;
    const double  hbarc = 197.329;


    if  ( n < 0)    nrerror( "Bad argument n in radial_harmonic_2d().");

    alpha = sqrt( mu*omega)/hbarc;
    if  ( alpha == 0.0 )    return  0.0;


    if  ( m < 0)    mm = -m;
    else    mm = m;


/*    Calculate energy.    */

    *E_nm = omega*( 2.0*n + 1.0*mm + 1.0);


/*    Calculate normalization.    */

    for  ( factor1 = 1, i = 1; i <= n; i++)    factor1 *= i;
    for  ( factor2 = 1, i = 1; i <= ( n + mm); i++)    factor2 *= i;
    norm = sqrt( alpha*alpha*( 2.0*factor1/factor2)/sqrt( PI));


/*    Calculate radial wavefunction.    */

    x = alpha*p;

    P_nm = norm;
    for  ( i = 1; i <= mm; i++)    P_nm *= x;
    P_nm *= exp( -x*x/2.0);
    P_nm *= laguerre_poly( n, 1.0*mm, x*x);

    return  P_nm;
}





/*******************************************************************************

        RADIAL_HARMONIC_3D

    This routine returns the 2-dimensional harmonic oscillator wave 
    function of order order n: H_n.

*******************************************************************************/

double  radial_harmonic_3d( double mu, double omega, int n, int l, 
                            double r, double *E_nl)

{
    int  i, factor0, factor1, factor2;
    double  alpha, x, norm, R_nl;
    const double  hbarc = 197.329;


    if  ( n < 1)    nrerror( "Bad argument n in radial_harmonic_3d().");
    if  ( l < 0)    nrerror( "Bad argument l in radial_harmonic_3d().");

    alpha = sqrt( mu*omega)/hbarc;
    if  ( alpha == 0.0 )    return  0.0;


/*    Calculate energy.    */

    *E_nl = omega*( 2.0*n + 1.0*l - 0.5);


/*    Calculate normalization.    */

    for  ( factor0 = 1, i = 1; i <= ( n + l + 1); i++)    factor0 *= 2;
    for  ( factor1 = 1, i = 1; i <= ( n - 1); i++)    factor1 *= i;
    for  ( factor2 = 1, i = 1; i <= ( 2*n + 2*l - 1); i += 2)    factor2 *= i;
    norm = sqrt( alpha*alpha*alpha*( 1.0*factor0*factor1/factor2)/sqrt( PI));


/*    Calculate radial wavefunction.    */

    x = alpha*r;

    R_nl = norm;
    for  ( i = 1; i <= l; i++)    R_nl *= x;
    R_nl *= exp( -x*x/2.0);
    R_nl *= laguerre_poly( n - 1, 1.0*l + 0.5, x*x);

    return  R_nl;
}

