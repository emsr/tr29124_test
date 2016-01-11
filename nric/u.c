
#include <stdio.h>
#include <math.h>


#include "nric.h"


main() {

    double ans, a, b, c;

    a = 2.0;
    b = -5.0;
    c = 7.0;


    printf( "\n\n  Test of utility routines." );

    printf( "\n\n" );

    ans = dsqr( a );
    printf( "\n  a = %lf  ans = dsqr( a ) = %lf", a, ans );

    ans = dcub( a );
    printf( "\n  a = %lf  ans = dcub( a ) = %lf", a, ans );

    ans = dpow4( a );
    printf( "\n  a = %lf  ans = dpow4( a ) = %lf", a, ans );

    ans = dpow5( a );
    printf( "\n  a = %lf  ans = dpow5( a ) = %lf", a, ans );

    ans = dsgn( a );
    printf( "\n  a = %lf  ans = dsgn( a ) = %lf", a, ans );

    ans = dsgn( b );
    printf( "\n  b = %lf  ans = dsgn( b ) = %lf", b, ans );

    ans = dsign( a, b );
    printf( "\n  a = %lf,  b = %lf,  ans = sign( a, b ) = %lf", a, b, ans );

    ans = dpythag( a, b );
    printf( "\n  a = %lf,  b = %lf,  ans = pythag( a, b ) = %lf", a, b, ans );

    ans = dmin( a, b );
    printf( "\n  a = %lf,  b = %lf,  ans = dmin( a, b ) = %lf", a, b, ans );

    ans = dmax( a, b );
    printf( "\n  a = %lf,  b = %lf,  ans = dmax( a, b ) = %lf", a, b, ans );

    ans = dmin3( a, b, c );
    printf( "\n  a = %lf,  b = %lf,  c = %lf,  ans = dmin3( a, b, c ) = %lf", a, b, c, ans );

    ans = dmax3( a, b, c );
    printf( "\n  a = %lf,  b = %lf,  c = %lf,  ans = dmax3( a, b, c ) = %lf", a, b, c, ans );



    printf( "\n\n" );
}


