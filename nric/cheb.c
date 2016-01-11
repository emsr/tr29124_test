
#include <stdio.h>
#include <math.h>

#include "nric.h"


double one( double x ) { return 1.0; }
double one_der( double x ) { return 0.0; }
double one_int( double x ) { return x; }

double ex( double x ) { return x; }
double ex_der( double x ) { return 1.0; }
double ex_int( double x ) { return x*x/2.0; }

double foo(double x) { return (1.0 - x)*exp(-x/2.0); }
double foo_der(double x) { return -0.5*(3.0 - x)*exp(-x/2.0); }
double foo_int(double x) { return 2.0*(1.0 + x)*exp(-x/2.0) - 2.0; }


main() {

    int i, n, m5, ma, mf;
    double a, b, h, x, f, err;
    double *c, *cint, *cder;

    a = 0.0;
    b = 5.0;
    h = (b-a)/10.0;
    n = 40;
    m5 = 5;
    ma = 10;
    mf = 15;

    c = dvector( 0, n-1 );
    cder = dvector( 0, n-1 );
    cint = dvector( 0, n-1 );


    printf( "\n\n\n  Test Chebyshev fitting.\n" );



    chebyshev_fit( a, b, c, n, one );
    chebyshev_deriv( a, b, c, cder, mf );
    chebyshev_integ( a, b, c, cint, mf );

    printf( "\n\n  Fit of f(x) = 1 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, n );
    for ( i = 0; i < n/2; ++i ) printf( "\n  c[%2d] = %16.12f                c[%2d] = %16.12f", i, c[i], i+n/2, c[i+n/2] );

    printf( "\n\n  Evaluation of f(x) = 1 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, m5 );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, m5, x );
        err = f - one(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of f(x) = 1 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, ma, x );
        err = f - one(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of derivative of f(x) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cder, ma, x );
        err = f - one_der(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of integral of f(x) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cint, ma, x );
        err = f - one_int(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }



    chebyshev_fit( a, b, c, n, ex );
    chebyshev_deriv( a, b, c, cder, mf );
    chebyshev_integ( a, b, c, cint, mf );

    printf( "\n\n  Fit of f(x) = x from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, n );
    for ( i = 0; i < n/2; ++i ) printf( "\n  c[%2d] = %16.12f                c[%2d] = %16.12f", i, c[i], i+n/2, c[i+n/2] );

    printf( "\n\n  Evaluation of f(x) = x from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, m5 );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, m5, x );
        err = f - ex(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of f(x) = x from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, ma, x );
        err = f - ex(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of derivative of f(x) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cder, ma, x );
        err = f - ex_der(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of integral of f(x) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cint, ma, x );
        err = f - ex_int(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }



    chebyshev_fit( a, b, c, n, foo );
    chebyshev_deriv( a, b, c, cder, mf );
    chebyshev_integ( a, b, c, cint, n );

    printf( "\n\n  Fit of f(x) = (1 - x)exp(-x/2) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, n );
    for ( i = 0; i < n/2; ++i ) printf( "\n  c[%2d] = %16.12f                c[%2d] = %16.12f", i, c[i], i+n/2, c[i+n/2] );

    printf( "\n\n  Fit of f'(x) = -(3 - x)exp(-x/2)/2 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, n );
    for ( i = 0; i < n/2; ++i ) printf( "\n  cder[%2d] = %16.12f                cder[%2d] = %16.12f", i, cder[i], i+n/2, cder[i+n/2] );

    printf( "\n\n  Fit of F(x) = 2(1 + x)exp(-x/2) - 2 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, n );
    for ( i = 0; i < n/2; ++i ) printf( "\n  cint[%2d] = %16.12f                cint[%2d] = %16.12f", i, cint[i], i+n/2, cint[i+n/2] );

    printf( "\n\n  Evaluation of f(x) = (1 - x)exp(-x/2) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, m5 );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, m5, x );
        err = f - foo(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of f(x) = (1 - x)exp(-x/2) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, ma );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, ma, x );
        err = f - foo(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of f(x) = (1 - x)exp(-x/2) from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, mf );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, c, mf, x );
        err = f - foo(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation of f'(x) = -(3 - x)exp(-x/2)/2 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, mf );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cder, mf, x );
        err = f - foo_der(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }

    printf( "\n\n  Evaluation F(x) = 2(1 + x)exp(-x/2) - 2 from a = %3.1f to b = %3.1f with %d coefficients . . .", a, b, mf );
    x = a;
    while ( x <= b ) {
        f = chebyshev_eval( a, b, cint, mf, x );
        err = f - foo_int(x);
        printf( "\n  f(%3.1f) = %16.12f    error = %16.12f", x, f, err );
        x += h;
    }




    printf( "\n\n" );


    free_dvector( cint, 0, n-1 );
    free_dvector( cder, 0, n-1 );
    free_dvector( c, 0, n-1 );
}



