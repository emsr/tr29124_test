

#include <math.h>


#include "nric.h"


/*
 *    Given a function funk and an initial guessd range x1 to x2, the routine
 *    expands the range geometrically until  a root is bracketed by the
 *    returned values x1 and x2 (in which case bracket returns 1)
 *    or until the range becomes unacceptably large ( in which case zbrak
 *    returns 0).  Success is guaranteed for a function which has opposite signs for
 *    sufficiently large and small arguments.
 */
int root_bracket( double (*funk)(double), double *x1, double *x2 ) {

    int i;
    double f1, f2;

    const double FACTOR = GOLD;
    const int NTRY = 50;

    if ( *x1 == *x2 ) nrerror( "Bad initial range in bracket." );
    f1 = (*funk)(*x1);
    f2 = (*funk)(*x2);
    for ( i = 1; i <= NTRY; ++i ) {
        if ( f1*f2 < 0.0 ) return 1;
        if ( fabs(f1) < fabs(f2) ) f1 = (*funk)( *x1 += FACTOR*( *x1 - *x2 ) );
        else f2 = (*funk)( *x2 += FACTOR*( *x2 - *x1 ) );
    }
    return 0;
}


/*
 *    Given a function funk defined on an interval x1 to x2, the routine
 *    subdivides the interval into n equally spaced segments, and searches
 *    for zero crossings of the function.  nb is the maximum number of roots
 *    sought,  and is reset to the number of bracketing pairs
 *    xb1[1..nb], xb2[1..nb] that are found.
 */
void root_brackets( double (*funk)(double),
               double x1, double x2, int n,
               double *xb1, double *xb2, int *nb ) {

    int nbb, i;
    double x, fp, fc, dx;

    nbb = 0;
    dx = (x2 - x1)/n;
    fp = (*funk)(x = x1);
    for ( i = 1; i <= n; ++i ) {
        fc = (*funk)(x += dx);
        if ( fc*fp <= 0.0 ) {
            xb1[++nbb] = x - dx;
            xb2[nbb] = x;
            if ( *nb == nbb ) return;
        }
        fp = fc;
    }
    *nb = nbb;

    return;
}



/*
 *    Using bisection, find the root of a function funk known to lie between x1 and x2.
 *    The root (which is returned) is refined until its accuracy is +/-eps.
 */
double root_bisect( double (*funk)(double), double x1, double x2, double eps ) {

    int i;
    double dx, f, fmid, xmid, x;

    const int IMAX = 55;

    f = (*funk)(x1);
    fmid = (*funk)(x2);
    if ( f*fmid >= 0.0 ) nrerror( "Root must be bracketed for bisection in root_bisect." );
    /*
     *    Orient search so that f > 0.0 lies at x + dx.
     */
    x = f < 0.0 ? (dx = x2-x1, x1) : (dx = x1-x2, x2);
    for ( i = 1; i <= IMAX; ++i ) {
        fmid = (*funk)( xmid = x + (dx *= 0.5) );
        if ( fmid < 0.0 ) x = xmid;
        if ( fabs(dx) < eps || fmid == 0.0 ) return x;
    }
    nrerror( "Too many bisections in root_bisect." );

    return  0.0;
}



/*
 *    Using the secant method, find the root of a function funk thought to
 *    lie between x1 and x2.  The root, returned as secant, is refined
 *    until its accuracy is +/- eps.
 */
double root_secant( double (*funk)(double), double x1, double x2, double eps ) {

    int i;
    double fl, f, dx, swap, xl, x;

    const double IMAX = 40;

    fl = (*funk)(x1);
    f = (*funk)(x2);
    if ( fabs(fl) < fabs(f) ) {
        x = x1;
        xl = x2;
        dswap( &fl, &f );
    }  else  { 
        xl = x;
        x = x2;
    }
    for ( i = 1; i <= IMAX; ++i ) {
        dx = ( xl - x )*f/( f - fl );
        xl = x;
        fl = f;
        x += dx;
        f = (*funk)(x);
        if ( fabs(dx) < eps || f == 0.0 ) return x;
    }
    nrerror( "Maximum number of iterations exceeded in root_secant." );

    return  0.0;
}



/*
 *    Using the false position method, find the root of a function funk known
 *    to lie between x1 and x2.  The root, returned as root_false_pos, is refined
 *    until its accuracy is +/- eps.
 */
double root_false_pos( double (*funk)(double), double x1, double x2, double eps ) {

    int i;
    double fl, fh, xl, xh, dx, swap, del, f, x;

    const double IMAX = 40;

    fl = (*funk)(x1);
    fh = (*funk)(x2);
    if ( fl*fh > 0.0 ) nrerror( "Root must be bracketed in root_false_pos." );
    if ( fl < 0.0 ) {
        xl = x1;
        xh = x2;
    }  else  { 
        xl = x2;
        xh = x1;
        dswap( &fl, &fh );
    }
    dx = xh - xl;
    for ( i = 1; i <= IMAX; ++i ) {
        x = xl + dx*fl/( fl - fh );
        f = (*funk)(x);
        if ( f < 0.0 ) {
            del = xl - x;
            xl = x;
            fl = f;
        }  else  {
            del = xh - x;
            xh = x;
            fh = f;
        }
        dx = xh - xl;
        if ( fabs(del) < eps || f == 0.0 )  return x;
    }
    nrerror( "Maximum number of iterations exceeded in root_false_pos." );

    return  0.0;
}


/*
 *    Using Ridder's method, find the root of a function funk known
 *    to lie between x1 and x2.  The root, returned as root_brent, is refined
 *    until its accuracy is +/- eps.
 */
double root_ridder( double (*funk)(double), double x1, double x2, double eps ) {

    int i;
    double ans, fh = (*funk)(x2), fl = (*funk)(x1), fm, fnew, s, xh = x2, xl = x1, xm, xnew;

    const int IMAX = 100;
    const double UNUSED = -1.0e30; /* an exceedingly unlikely answer. */

    if ( fl*fh < 0.0 ) {
        ans = UNUSED;
        for ( i = 1; i <= IMAX; ++i ) {
            xm = 0.5*(xl + xh);
            fm = (*funk)(xm);
            s = sqrt(fm*fm - fl*fh);
            if ( s == 0.0 ) return ans;
            xnew = xm + (xm - xl)*( fl >= fh ? 1.0 : -1.0 )*fm/s;
            if ( fabs(xnew - ans) < eps ) return ans;
            fnew = (*funk)(ans = xnew);
            if ( fnew == 0.0 ) return ans;
            if ( dsign( fm, fnew) != fm) {
                xl = xm;
                fl = fm;
                xh = xnew;
                fh = fnew;
            } else if ( dsign(fl, fnew) != fl ) {
                xh = xnew;
                fh = fnew;
            } else if ( dsign(fh, fnew) != fh ) {
                xl = xnew;
                fl = fnew;
            } else nrerror( "Some major malfunction in root_ridder." );
            if ( fabs(xh - xl) < eps ) return ans;
        }
        nrerror( "Maximum number of iterations exceeded in root_ridder." );
    } else {
        if ( fl == 0.0 ) return x1;
        if ( fh == 0.0 ) return x2;
        nrerror( "Root must be bracketed in root_ridder." );
    }

    return  0.0;
}



/*
 *    Using Brent's method, find the root of a function funk known to lie between x1 and x2.
 *    The root, returned as brent, will be refined until it's accuracy is eps.
 */
double root_brent( double (*funk)( double), double x1, double x2, double eps ) {

    int iter;
    double a = x1, b = x2, c=x2, d, e, min1, min2;
    double fa = (*funk)(a), fb = (*funk)(b), fc, p, q, r, s, tol1, xm;

    const int ITMAX = 100;
    const double EPS = 1.0e-12;

    if ( fb*fa > 0.0 ) nrerror( "Root must be bracketed in root_brent." );
    fc = fb;
    for ( iter = 1; iter <= ITMAX; ++iter ) {

        if ( fb*fc > 0.0 ) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if ( fabs(fc) < fabs(fb) ) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        tol1 = 2.0*EPS*fabs(b) + 0.5*eps;
        xm = 0.5*( c - b );
        if ( fabs(xm) <= tol1 || fb == 0.0 ) return b;
        if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) ) {
            s = fb/fa;
            if ( a == c ) {
                p = 2.0*xm*s;
                q = 1.0 - s;
            } else {
                q = fa/fc;
                r = fb/fc;
                p = s*( 2.0*xm*q*( q - r ) - ( b - a )*( r - 1.0 ) );
                q = ( q - 1.0 )*( r - 1.0 )*( s - 1.0 );
            }
            if ( p > 0.0 ) q = -q;
            p = fabs(p);
            min1 = 3.0*xm*q - fabs( tol1*q );
            min2 = fabs( e*q );
            if ( 2.0*p < dmin( min1, min2 ) ) {
                e = d;
                d = p/q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if ( fabs(d) > tol1 ) b += d;
        else b += dsign( tol1, xm );
        fb = (*funk)(b);
    }
    nrerror( "Maximum number of iterations exceeded in root_brent." );

    return  0.0;
}



/*
 *    Using the Newton-Raphson method, find the root of a function known to lie in the interval
 *    x1 to x2.  The root will be refined until its accuracy is known within +/- eps.
 *    funcd is a user-supplied routine that provides both the function and the first derivative
 *    of the function at the point x.
 */
double root_newton( void (*funky)( double, double *, double * ), double x1, double x2, double eps ) {

    int i;
    double df, dx, f, x;

    const double IMAX = 40;

    x = 0.5*( x1 + x2 );
    for ( i = 1; i <= IMAX; ++i ) {
        (*funky)( x, &f, &df );
        dx = f/df;
        x -= dx;
        if ( (x1 - x)*(x - x2) < 0.0 ) nrerror( "Jumped out of brackets in root_newton." );
        if ( fabs(dx) < eps ) return x;
    }
    nrerror( "Maximum number of iterations in root_newton." );

    return  0.0;
}



/*
 *    Using a combination of Newton-Raphson and bisection, find the root of a function known to lie in the interval
 *    x1 to x2.  The root will be refined until its accuracy is known within +/- eps.
 *    funcd is a user-supplied routine that provides both the function and the first derivative
 *    of the function at the point x.
 */
double root_safe( void (*funky)( double, double *, double * ), double x1, double x2, double eps ) {

    int i;
    double df, dx, dxold, f, fh, fl, temp, xh, xl, x;

    const int IMAX = 100;

        (*funky)( x1, &fl, &df );
        (*funky)( x2, &fh, &df );
    if ( fl*fh < 0.0 ) nrerror( "Root must be bracketed in root_safe." );
    if ( fl == 0.0 ) return x1;
    if ( fh == 0.0 ) return x2;
    if ( fl < 0.0 ) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    x = 0.5*(x1 + x2);
    dxold = fabs(x2 - x1);
    dx = dxold;
    for ( i = 1; i <= IMAX; ++i ) {
        if ( ((x - xh)*df - f)*((x - xl)*df - f) >= 0.0 || fabs(2.0*f) > fabs(dxold*df) ) {
            dxold = dx;
            dx = 0.5*(xh - xl);
            x = xl + dx;
            if ( x = xl ) return x;
        } else {
            dxold = dx;
            dx = f/df;
            temp = x;
            x -= dx;
            if ( temp == x ) return x;
        }
        if ( fabs(dx) < eps ) return x;
        (*funky)( x,  &f, &df );
        if ( f < 0.0 ) xl = x;
        else xh = x;
    }
    nrerror( "Maximum number of iterations in root_safe." );

    return  0.0;
}



