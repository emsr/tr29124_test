

#include <math.h>

#include "nric.h"



/*
 *    Trapeziod rule integration.
 *
 *    This routine implements the nth stage of refinement of an extended trapezoid rule.
 *    With n = 1, the crudest estimate of the integral is returned.
 *    With successive calls with n = 2, 3, ... (in order) accuracy will be improved
 *    by adding 2^(n-2) interior points.
 */

double trapezoid( double (*funk)(double), double a, double b, int n ) {

    double x, sum, del;
    int j;
    static double s;
    static int lastn, it;

    if ( n <= 0 ) nrerror( "Non-positive order in trapezoid." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in trapezoid." );
    lastn = n;

    if ( n == 1 ) {
        s = 0.5*( b - a )*( (*funk)(a) + (*funk)(b) );
        it = 1;
    } else {
        del = (b - a)/it;
        x = a + 0.5*del;
        for ( sum = 0.0, j = 1; j <= it; ++j, x += del ) sum += (*funk)(x);
        s = 0.5*( s + ( b - a )*sum/it );
        it *= 2;
    }
    return s;
}



/*
 *    Modified midpoint integration.
 *
 *    This routine implements the nth stage of refinement of a modified midpoint integration.
 *    With n = 1, the crudest estimate of the integral is returned.
 *    With successive calls with n = 2, 3, ... (in order) accuracy will be improved
 *    by adding (2/3)^(n-1) interior points.
 */

double midpoint( double (*funk)(double), double a, double b, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return s = (b - a)*(*funk)(x);
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += (*funk)(x);
            x += ddel;
            sum += (*funk)(x);
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s + (b - a)*sum/tnm)/3.0;
    }
}


/*
 *    This routine is an exact replacement of midpoint except that the
 *    points are evenly spaced in 1/x rather than x.  This allows the
 *    lower limit aa to be as large and negative or the upper limit
 *    bb to be as large and positive as the computer allows but not both.
 *    aa and bb must have the same sign.
 */
double midpoint_inv( double (*funk)(double), double aa, double bb, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel, a, b;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    a = 1.0/aa;
    b = 1.0/bb;

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return s = (b - a)*(*funk)(1.0/x)/(x*x);
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += (*funk)(1.0/x)/(x*x);
            x += ddel;
            sum += (*funk)(1.0/x)/(x*x);
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s + (b - a)*sum/tnm)/3.0;
    }
}


/*
 *    This routine is an exact replacement of midpoint except that it allows
 *    for an inverse square root singularity at the upper limit bb.
 */
double midpoint_inv_sqrt_lower( double (*funk)(double), double aa, double bb, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel, a, b;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    a = 0.0;
    b = sqrt(bb - aa);

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return s = (b - a)*2.0*x*(*funk)(aa + x*x);
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += 2.0*x*(*funk)(aa + x*x);
            x += ddel;
            sum += 2.0*x*(*funk)(aa + x*x);
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s + (b - a)*sum/tnm)/3.0;
    }
}



/*
 *    This routine is an exact replacement of midpoint except that it allows
 *    for an inverse square root singularity at the upper limit bb.
 */
double midpoint_inv_sqrt_upper( double (*funk)(double), double aa, double bb, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel, a, b;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    a = 0.0;
    b = sqrt(bb - aa);

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return s = (b - a)*2.0*x*(*funk)(bb - x*x);
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += 2.0*x*(*funk)(bb - x*x);
            x += ddel;
            sum += 2.0*x*(*funk)(bb - x*x);
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s + (b - a)*sum/tnm)/3.0;
    }
}



/*
 *    This routine is an exact replacement of midpoint except that bb is assumed
 *    to be infinite (input value is not actually used).  It is assumed that the
 *    function funk decreases exponentially rapidly at infinity.
 */
double midpoint_exp( double (*funk)(double), double aa, double bb, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel, a, b;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    a = 0.0;
    b = exp(-aa);

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return (s = (b - a)*(*funk)(-log(x))/x);
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += (*funk)(-log(x))/x;
            x += ddel;
            sum += (*funk)(-log(x))/x;
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s + (b - a)*sum/tnm)/3.0;
    }
}



/*
 *    This routine is an exact replacement of midpoint except that aa is assumed
 *    to be zero (input value is not actually used).  It is assumed that the
 *    function funk has a logarithmic singularity at 0.
 */
double midpoint_log( double (*funk)(double), double aa, double bb, int n ) {

    int j, it;
    double x, tnm, sum, del, ddel, a, b;
    static int lastn;
    static double s;

    if ( n <= 0 ) nrerror( "Non-positive order in midpoint." );
    if ( n != 1 && n != lastn + 1 ) nrerror( "Order out of sequence in midpoint." );
    lastn = n;

    a = -log(bb);
    b = -log(0.01);

    if ( n == 1 ) {
        x = 0.5*(a + b);
        return (s = -(b - a)*x*(*funk)(exp(-x)));
    } else {
        for ( it = 1, j = 1; j < n-1; ++j ) it *= 3;
        tnm = it;
        /*
         *    The added points alternate in spacing between del and ddel.
         */
        del = (b - a)/(3.0*tnm);
        ddel = 2.0*del;
        x = a + 0.5*del;
        sum = 0.0;
        for ( j = 1; j <= it; ++j ) {
            sum += x*(*funk)(exp(-x));
            x += ddel;
            sum += x*(*funk)(exp(-x));
            x += del;
        }
        /*
         *    The new sum is combined with the old integral to give a refined integral.
         */
        return s = (s - (b - a)*sum/tnm)/3.0;
    }
}



/*
 *    Runs through n steps of the trapezoid rule integration
 *    of a function funk of one real variable from a to b.
 */
double dumb_trapezoid( double (*funk)(double), double a, double b, int n ) {

    int j;
    double s;
    const int JMAX = 20;

    if ( n <= 0 ) nrerror( "Non-positive order in dumb_trapezoid." );
    if ( n > JMAX ) nrerror( "Order too large in dumb_trapezoid." );

    for ( j = 1; j <= n; ++j ) { s = trapezoid( funk, a, b, j ); }
    return s;
}



/*
 *    Integrates a function of one real variable over the interval a to b
 *    using trapezoid rule integration.  Integration steps are taken until
 *    the difference between successive steps is less than eps.
 */
double quad_trapezoid( double (*funk)(double), double a, double b, double eps ) {

    int j;
    double s, olds;
    const int JMAX = 20;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in quad_trapezoid." );

    olds = -1.0e30;
    for ( j = 1; j <= JMAX; ++j ) {
        s = trapezoid( funk, a, b, j );
        if ( fabs(s - olds) < eps*fabs(olds) ) return s;
        if ( fabs(s) < eps && fabs(olds) < eps && j > 6 ) return s;
        olds = s;
    }
    nrerror( "Too many steps in routine quad_trapezoid." );
    return 0.0;
}


/*
 *    Runs through n steps of the Simpson rule integration
 *    of a function funk of one real variable from a to b.
 */
double dumb_simpson( double (*funk)(double), double a, double b, int n ) {

    int j;
    double s, st, ost;
    const int JMAX = 20;

    if ( n <= 0 ) nrerror( "Non-positive order in dumb_simpson." );
    if ( n > JMAX ) nrerror( "Order too large in dumb_simpson." );

    for ( j = 1; j <= n; ++j ) {
        st = trapezoid( funk, a, b, j );
        if ( j == 1 ) ost = st;
        s = ( 4.0*st - ost )/3.0;
        ost = st;
    }
    return s;
}



/*
 *    Integrates a function of one real variable over the interval a to b
 *    using Simpson rule integration.  Integration steps are taken until
 *    the difference between successive steps is less than eps.
 */
double quad_simpson( double (*funk)(double), double a, double b, double eps ) {

    int j;
    double s, olds, st, oldst;
    const int JMAX = 20;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in quad_simpson." );

    oldst = olds = -1.0e30;
    for ( j = 1; j <= JMAX; ++j ) {
        st = trapezoid( funk, a, b, j );
        s = ( 4.0*st - oldst )/3.0;
        if ( fabs(s - olds) < eps*fabs(olds) ) return s;
        if ( fabs(s) < eps && fabs(olds) < eps && j > 6 ) return s;
        olds = s;
        oldst = st;
    }
    nrerror( "Too many steps in routine quad_simpson." );
    return 0.0;
}


/*
 *    Runs through n steps of the Romberg integration
 *    of a function funk of one real variable from a to b.
 */
double dumb_romberg( double (*funk)(double), double a, double b, int n ) {

    int j;
    double ss, dss, *s, *h;
    const int K = 5;
    const int JMAX = 20;

    if ( n <= 0 ) nrerror( "Non-positive order in dumb_romberg." );
    if ( n > JMAX ) nrerror( "Order too large in dumb_romberg." );

    s = dvector( 0, JMAX );
    h = dvector( 0, JMAX+1 );

    h[1] = 1.0;
    for ( j = 1; j <= n; ++j ) {
        s[j] = trapezoid( funk, a, b, j );
        if ( j >= K )  {
            poly_interp( &h[j-K], &s[j-K], K, 0.0, &ss, &dss );
        }
        h[j+1] = 0.25*h[j];
    }

    free_dvector( h, 1, JMAX+1 );
    free_dvector( s, 1, JMAX );

    return ss;
}



/*
 *    Integrates a function of one real variable over the interval a to b
 *    using Romberg integration.  Integration steps are taken until the
 *    difference between successive steps is less than eps.
 */
double quad_romberg( double (*funk)(double), double a, double b, double eps ) {

    int j;
    double ss, dss, *s, *h;
    const int K = 5;
    const int JMAX = 20;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in quad_romberg." );

    s = dvector( 1, JMAX );
    h = dvector( 1, JMAX+1 );

    h[1] = 1.0;
    for ( j = 1; j <= JMAX; ++j ) {
        s[j] = trapezoid( funk, a, b, j );
        if ( j >= K )  {
            poly_interp( &h[j-K], &s[j-K], K, 0.0, &ss, &dss );
            if ( fabs(dss) < eps*fabs(ss) ) return ss;
            if ( fabs(dss) < eps && fabs(ss) < eps && j > 6 ) return ss;
        }
        h[j+1] = 0.25*h[j];
    }

    free_dvector( h, 1, JMAX+1 );
    free_dvector( s, 1, JMAX );

    nrerror( "Too many steps in routine quad_romberg." );
    return 0.0;
}




/*
 *    Romberg integration on an open interval.  Returns the integral of a function
 *    funk from a to b using any specified integration routine choose and Romberg's method.
 *    Normally, choose will be an open formula, not evaluating the function at the endpoints.
 *    It is assumed that choose triples the number of steps on each call, and that
 *    it's error series contains only even powers of the steps.  The routines
 *    midpoint, midpoint_inv, midpoint_inv_sqrt_lower, midpoint_inv_sqrt_upper, and
 *    midpoint_exp are possible choices for choose.
 */
double quad_romberg_open( double (*funk)(double), double a, double b, double eps,
                       double (*choose)( double (*)(double), double, double, int ) ) {

    int j;
    double ss, dss, *s, *h;
    const int K = 5;
    const int JMAX = 14;

    if ( eps <= 0.0 ) nrerror( "Error tolerance eps must be greater than 0 in quad_romberg_open." );

    s = dvector( 1, JMAX );
    h = dvector( 1, JMAX+1 );

    h[1] = 1.0;
    for ( j = 1; j <= JMAX; ++j ) {
        s[j] = (*choose)( funk, a, b, j );
        if ( j >= K )  {
            poly_interp( &h[j-K], &s[j-K], K, 0.0, &ss, &dss );
            if ( fabs(dss) < eps*fabs(ss) ) return ss;
            if ( fabs(dss) < eps && fabs(ss) < eps && j > K+1 ) return ss;
        }
        h[j+1] = h[j]/9.0;
    }

    free_dvector( h, 1, JMAX+1 );
    free_dvector( s, 1, JMAX );

    nrerror( "Too many steps in routine quad_romberg_open." );
    return 0.0;
}


