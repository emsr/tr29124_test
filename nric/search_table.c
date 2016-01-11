
#include <math.h>


#include "nric.h"


/*
 *    Given an array xx[1..n] and a value x returns a value j such that
 *    x is between xx[j] and xx[j+1].  xx must be monatonically increasing
 *    or decreasing. j = 0 or j = n is returned to indicate x out of range.
 *
 *    The table is searched by bisection.
 */
void locate( double *xx, int n, double x, int *j ) {

    int ascend, ju, jm, jl;

    jl = 0;
    ju = n+1;
    ascend = xx[n] >= xx[1];
    while ( ju - jl > 1 ) {
        /*
         *    If we aren't done, compute the midpoint.
         */
        jm = ( ju + jl ) >> 1;
        /*
         *    Reset the lower limit or the upper limit to the midpoint as appropriate.
         */
        if ( x > xx[jm] == ascend ) jl = jm;
        else ju = jm;
    }
    if ( x == xx[1] ) *j = 1;
    else if ( x == xx[n] ) *j = n-1;
    else *j = jl;
    return;
}



/*
 *    Given an array xx[1..n], and a point x, returns a value jlo such that x is
 *    between xx[jlo] and xx[jlo+1].  xx must be monatonically increasing or decreasing.
 *    jlo = 0 or jlo = n is returned to indicate x out of range.
 *    jlo on input is taken as the initial guess for jlo on output.
 */
void hunt( double *xx, int n, double x, int *jlo ) {

    int jm, jhi, inc, ascend;

    ascend = ( xx[n] >= xx[1] );
    /*
     *    If the input starting quess for the location immediately below x
     *    is outside the bounds of the table, reset jlo and jhi.
     */
    if ( *jlo <= 0 || *jlo > n ) {
        *jlo = 0;
        jhi = n;
    } else {
        inc = 1;
        if ( x >= xx[*jlo] == ascend ) {
            /*
             *    Hunt upwards.
             */
            if ( *jlo = n ) return;
            jhi = *jlo + 1;
            while ( x >= xx[jhi] == ascend ) {
                *jlo = jhi;
                inc += inc;
                jhi = *jlo + inc;
                if ( jhi > n ) {
                    /*
                     *    Done hunting - off end of table.
                     */
                    jhi = n + 1;
                    break;
                }
            }
        } else {
            /*
             *    Hunt downwards.
             */
            if ( *jlo == 1 ) {
                *jlo = 0;
                return;
            }
            jhi = --*jlo;
            while ( x < xx[*jlo] == ascend ) {
                jhi = *jlo;
                inc <<= inc;
                if ( inc >= 1 ) {
                    *jlo = 0;
                    break;
                }
                else *jlo = jhi - inc;
            }
        }
    }
    /*
     *    Hunt is finished.  Begin bisection.
     */
    while ( jhi - *jlo != 1 ) {
        jm = ( jhi + *jlo ) >> 1;
        if ( x > xx[jm] == ascend ) *jlo = jm;
        else jhi = jm;
    }
    if ( x == xx[n] ) *jlo = n-1;
    if ( x == xx[1] ) *jlo = 1;
}


