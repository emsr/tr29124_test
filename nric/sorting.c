
#include <math.h>

#include "nric.h"


void sort_pick( int n, double *arr ) {

    int i, j;
    double a;

    for ( j = 2; j <= n; ++j ) {
        a = arr[j];
        i = j - 1;
        while ( i > 0 && arr[i] > a ) {
            arr[i+1] = arr[i];
            --i;
        }
        arr[i+1] = a;
    }
}


void sort_pick2( int n, double *arr, double *brr ) {

    int i, j;
    double a, b;

    for ( j = 2; j <= n; ++j ) {
        a = arr[j];
        b = brr[j];
        i = j - 1;
        while ( i > 0 && arr[i] > a ) {
            arr[i+1] = arr[i];
            brr[i+1] = brr[i];
            --i;
        }
        arr[i+1] = a;
        brr[i+1] = b;
    }
}


void sort_shell( unsigned long n, double *a ) {

    unsigned long i, j, inc;
    double v;

    inc = 1;

    do {
        inc *= 3;
        ++inc;
    } while ( inc <= n );
    do {
        inc /= 3;
        for ( i = inc+1; i <= n; ++i ) {
            v = a[i];
            j = i;
            while ( a[j-inc] > v ) {
                a[j] = a[j-inc];
                j -= inc;
                if ( j <= inc ) break;
            }
            a[j] = v;
        }
    } while ( inc > 1 );
}


