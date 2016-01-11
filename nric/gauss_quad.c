

#include <math.h>


#include "nric.h"

#define EPS 1.0e-15
#define MAXIT 40
#define PIM4 0.7511255444649425


/*
 *    This sums up the the wights times function evaluations for gauss quadrature formulas.
 *
 *    Warning ! ! !
 *
 *    This assumes that the abscissas have already been scaled from -1 < y < +1
 *    to a < x < b and that the weights have been scaled by (b - a)/2r.
 *    This cannot be used with the other quadratures that use standard unscaled
 *    abscissas and weights.
 */
double quad_gauss_legendre( double (*funk)(double), double *x, double *w, int n ) {

    int i;
    double sum;

    for ( sum = 0.0, i = 1; i <= n; ++i ) sum += w[i]*funk(x[i]);
    return sum;
}


/*
 *    This will do the shifting of abscissas and weights
 *    from  -1 < y < +1  to  a < x < b  and
 *    from       w       to  (b - a)w/2.
 */
double quad_gauss( double (*funk)(double), double a, double b, double *x, double *w, int n ) {

    int i;
    double sum, bpa, bma;

    bpa = 0.5*(b + a);
    bma = 0.5*(b - a);

    for ( sum = 0.0, i = 1; i <= n; ++i ) sum += w[i]*funk(bpa + bma*x[i]);
    return bma*sum;
}



/*****************************************************************************************

        gauss_legendre

    This routine is taken from Numerical Recipies in C by Press, Flannery, Teukolsky, 
    and Vetterling.  It calculates wieghts and grid points for gaussian quadrature 
    integration.

*****************************************************************************************/

void gauss_legendre( int n, double a, double b, double *x, double *w) {

    int m, k = 0;
    double bpa, bma;
    int i, j;
    double z1, z;
    double p2, p1, p, pp;

    m = ( n + 1)/2;
    bpa = ( b + a)/2.0;
    bma = ( b - a)/2.0;

    for  ( i = 1; i <= m; i++)    {
        z = cos( PI*( i - 1.0/4.0)/( n + 1.0/2.0));    /*    Clever approximation of root.    */
        k = 0;
        do    {
            /*
             *    Compute p, p1, and p2 the Legendre polynomials of order n, n-1, n-2 respectively
             *    by iterating through the recursion relation for the Legendre polynomials.
             *    Compute pp the derivative of the Legendre polynomial of order n.
             */
            for  ( p = 1.0, p1 = 0.0, j = 1; j <= n; j++)    {
                p2 = p1;
                p1 = p;
                p = ( ( 2.0*j - 1.0)*z*p1 - ( j - 1.0)*p2)/j;  /*  Recursion relation for legendre polynomials.  */
            }
            pp = n*( z*p - p1)/( z*z - 1.0);  /*  Recursion relation for derivatives of legendre polynomials.  */
            z1 = z;
            z = z1 - p/pp;  /*  Converge on root by Newton's method.  */
            k = k + 1;
        } while  ( fabs( z - z1) > EPS);

        x[i] = bpa - bma*z;
        x[n+1-i] = bpa + bma*z;
        w[i] = 2.0*bma/( ( 1.0 - z*z)*pp*pp);
        w[n+1-i] = w[i];
    }
}



void gauss_laguerre( double *x, double *w, int n, double alpha ) {

    int i, its, j;
    double ai;
    double p1, p2, p3, pp, z, z1;

    for ( i = 1; i <= n; ++i ) {
        if ( i == 1 ) z = (1.0 + alpha)*(3.0 + 0.92*alpha)/(1.0 + 2.4*n + 1.8*alpha);
        else if ( i == 2 ) z += (15.0 + 6.25*alpha)/(1.0 + 2.5*n + 0.9*alpha);
        else {
            ai = i - 2;
            z += ((1.0 + 2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0 + 3.5*ai))*(z - x[i-2])/(1.0 + 0.3*alpha);
        }
        for ( its = 1; its <= MAXIT; ++its ) {
            p1 = 1.0;
            p2 = 0.0;
            for ( j = 1; j <= n; ++j ) {
                p3 = p2;
                p2 = p1;
                p1 = ((2*j - 1 + alpha - z)*p2 - (j - 1 + alpha)*p3)/j;
            }
            pp = (n*p1 - (n + alpha)*p2)/z;
            z1 = z;
            z = z1 - p1/pp;
            if ( fabs(z - z1) <= 100*EPS ) break;
        }
        if ( its > MAXIT ) nrerror( "Too many iterations in gauss_laguerre." );
        x[i] = z;
        w[i] = -exp( ln_gamma(alpha + n) - ln_gamma(1.0*n) )/(pp*n*p2);
    }
}



void gauss_hermite( double *x, double *w, int n ) {

    int i, its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n+1)/2;
    for ( i = 1; i <= m; ++i ) {
        if ( i == 1 ) z = sqrt( 2.0*n + 1.0 ) - 1.85575*pow( 2.0*n + 1.0, -0.166667 );
        else if ( i == 2 ) z -= 1.14*pow( 1.0*n, 0.426)/z;
        else if ( i == 3 ) z = 1.86*z - 0.86*x[1];
        else if ( i == 4 ) z = 1.91*z - 0.91*x[2];
        else z = 2.0*z - x[i-2];
        for ( its = 1; its <= MAXIT; ++its ) {
            p1 = PIM4;
            p2 = 0.0;
            for ( j = 1; j <= n; ++j ) {
                p3 = p2;
                p2 = p1;
                p1 = z*sqrt(2.0/j)*p2 - sqrt(1.0*(j-1)/j)*p3;
            }
            pp = sqrt(2.0*n)*p2;
            z1 = z;
            z = z1 - p1/pp;
            if ( fabs(z - z1) <= EPS ) break;
        }
        if ( its > MAXIT ) nrerror( "Too many iterations in gauss_hermite." );
        x[i] = z;
        x[n+1-i] = -z;
        w[i] = 2.0/(pp*pp);
        w[n+1-i] = w[i];
    }
}



void gauss_jacobi( double *x, double *w, int n, double alpha, double beta ) {

    int i, its, j;
    double alphabeta, an, bn, r1, r2, r3;
    double a, b, c, p1, p2, p3, pp, temp, z, z1;

    for ( i = 1; i <= n; ++i ) {
        if ( i == 1 ) {
            an = alpha/n;
            bn = beta/n;
            r1 = (1.0 + alpha)*(2.78/(4.0 + n*n) + 0.768*an/n);
            r2 = 1.0 + 1.48*an + 0.96*bn + 0.452*an*an + 0.83*an*bn;
            z = 1.0 - r1/r2;
        } else if ( i == 2 ) {
            r1 = (4.1 + alpha)/((1.0 + alpha)*(1.0 + 0.156*alpha));
            r2 = 1.0 + 0.06*(n - 8.0)*(1.0 + 0.12*alpha)/n;
            r3 = 1.0 + 0.012*beta*(1.0 + 0.25*fabs(alpha))/n;
            z -= (1.0 - z)*r1*r2*r3;
        } else if ( i == 3 ) {
            r1 = (1.67 + 0.28*alpha)/(1.0 + 0.37*alpha);
            r2 = 1.0 + 0.22*(n - 8.0)/n;
            r3 = 1.0 + 8.0*beta/((6.28 + beta)*n*n);
            z -= (x[1] - z)*r1*r2*r3;
        } else if ( i == n-1 ) {
            r1 = (1.0 + 0.235*beta)/(0.766 + 0.119*beta);
            r2 = 1.0/(1.0 + 0.639*(n - 4.0)/(1.0 + 0.71*(n - 4.0)));
            r3 = 1.0/(1.0 + 20.0*alpha/((7.5 + alpha)*n*n));
            z += (z - x[n-3])*r1*r2*r3;
        } else if ( i == n ) {
            r1 = (1.0 + 0.37*beta)/(1.67 + 0.28*beta);
            r2 = 1.0/(1.0 + 0.22*(n - 8.0)/n);
            r3 = 1.0/(1.0 + 8.0*alpha/((6.28 + alpha)*n*n));
            z += (z - x[n-2])*r1*r2*r3;
        } else  {
            z = 3.0*x[i-1] - 3.0*x[i-2] + x[i-3];
        }
        alphabeta = alpha + beta;
        for ( its = 1; its <= MAXIT; ++its ) {
            temp = 2.0 + alphabeta;
            p1 = (alpha - beta + temp*z)/2.0;
            p2 = 1.0;
            for ( j = 2; j <= n; ++j ) {
                p3 = p2;
                p2 = p1;
                temp = 2.0*j + alphabeta;
                a = 2.0*j*(j + alphabeta)*(temp - 2.0);
                b = (temp - 1.0)*(alpha*alpha - beta*beta + temp*(temp - 2.0)*z );
                c = 2.0*(j - 1 + alpha)*(j - 1 + beta)*temp; 
                p1 = (b*p2 - c*p3)/a;
            }
            pp = (n*(alpha - beta - temp*z)*p1 + 2.0*(n + alpha)*(n + beta)*p2)/(temp*(1.0 - z*z));
            z1 = z;
            z = z1 - p1/pp;
            if ( fabs(z - z1) <= EPS ) break;
        }
        if ( its > MAXIT ) nrerror( "Too many iterations in gauss_jacobi." );
        x[i] = z;
        w[i] = exp( ln_gamma(alpha + n)
                  + ln_gamma(beta + n)
                  - ln_gamma(n + 1.0)
                  - ln_gamma(n + alphabeta + 1.0))*temp*pow( 2.0, alphabeta)/(pp*p2);
    }
}


/*
 *    Generates the Gauss-Chebyshev abscissas and the weight (a constant).
 *    This gives the same results as gauss_jacobi for alpha = beta = -0.5.
 */
void gauss_chebyshev( double *x, double *w, int n ) {

    int i, m;

    m = (n + 1)/2;
    for ( i = 1; i <= m; ++i ) {
        x[n+1-i] = -(x[i] = cos(PI*(i-0.5)/n));
        w[n+1-i] = w[i] = PI/n;
    }
}


double gauss_crap( double (*funk)(double), double a, double b, int n ) {

    int j, jmax;
    double y, sum;
    static int nn, oldn;
    static double bpa, bma, oldsum;

    if ( n <= 0 ) nrerror( "Non-positive order in gauss_crap." );
    if ( n != 1 && n != oldn+1 ) nrerror( "Order out of sequence in gauss_crap." );

    if ( n == 1 ) {
        bpa = 0.5*(b + a);
        bma = 0.5*(b - a);
        nn = 3;
        y = cos(PI/6.0);
        oldn = n;
        return oldsum = PI*bma*((*funk)(bpa + y*bma) + (*funk)(bpa) + (*funk)(bpa - y*bma))/nn;
    } else {
        nn *= 3;
        y = cos(PI*0.5/nn);
        sum = (*funk)(bpa + y*bma) + (*funk)(bpa - y*bma);
        jmax = 1 + (nn - 4 - n)/6;
        for ( j = 1; j <= jmax; ++j ) {
            y = cos(PI*(3.0*j - 0.5)/nn);
            sum = (*funk)(bpa + y*bma) + (*funk)(bpa - y*bma);
            y = cos(PI*(3.0*j + 0.5)/nn);
            sum = (*funk)(bpa + y*bma) + (*funk)(bpa - y*bma);
        }
        oldn = n;
        return oldsum = oldsum/3.0 + PI*bma*sum/nn;
    }
}


double dumb_gauss_crap( double (*funk)(double), double a, double b, int n ) {

    int j;
    double s, olds;
    const int JMAX = 12;

    if ( n <= 0 ) nrerror( "Non-positive order in dumb_gauss_crap." );
    if ( n > JMAX ) nrerror( "Order too large in dumb_gauss_crap." );

    for ( j = 1; j <= n; ++j ) s = gauss_crap( funk, a, b, j );
    return s;
}


double quad_gauss_crap( double (*funk)(double), double a, double b, double eps ) {

    int j;
    double s, olds;
    const int JMAX = 12;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in quad_gauss_crap." );

    olds = -1.0e30;
    for ( j = 1; j <= JMAX; ++j ) {
        s = gauss_crap( funk, a, b, j );
        if ( fabs(s - olds) < eps*fabs(olds) ) return s;
        if ( fabs(s) < eps && fabs(olds) < eps && j > 6 ) return s;
        olds = s;
    }
    nrerror( "Too many steps in routine quad_gauss_crap." );
    return 0.0;
}


