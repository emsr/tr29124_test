

#include <math.h>

#include "nric.h"

#define MAXIT 100
#define EPS 1.0e-10
#define FPMIN 1.0e-30


double exp_int( int n, double x ) {

    int i, ii, nm1;
    double a, b, c, d, del, fact, h, psi, ans;

    nm1 = n-1;

    if ( n < 0 || x < 0.0 || ( x == 0.0 && ( n == 0 || n == 1 ) ) ) nrerror( "Bad arguments in exp_int." );

    if ( n == 0 ) return ans = exp(-x)/x;

    if ( x == 0.0 ) return ans = 1.0/nm1;

    if ( x > 1.0 ) {
        b = x + n;
        c = 1.0/FPMIN;
        d = 1.0/b;
        h = d;
        for ( i = 1; i <= MAXIT; ++i ) {
            a = -i*(nm1 + i);
            b += 2.0;
            d = 1.0/(a*d + b);
            c = b + a/c;
            del = c*d;
            h *= del;
            if ( fabs(del - 1.0) < EPS ) {
                ans = h*exp(-x);
                return ans;
            }
        }
        nrerror( "Continued fraction failed in exp_int." );
    } else {
        ans = (nm1 != 0 ? 1.0/nm1 : -log(x) - GAMMA );
        fact = 1.0;
        for ( i = 1; i <= MAXIT; ++i ) {
            fact *= -x/i;
            if ( i != nm1 ) del = -fact/(i - nm1);
            else {
                psi = -GAMMA;
                for ( ii = 1; ii <= nm1; ++ii ) psi += 1.0/ii;
                del = fact*(psi - log(x)); 
            }
            ans += del;
            if ( fabs(del) < EPS*fabs(ans) ) return ans;
        }
        nrerror( "Series failed in exp_int." );
    }
    return ans;
}


double ei( double x ) {

    int k;
    double fact, prev, sum, term;

    if ( x <= 0.0 ) nrerror( "Bad argument in ei." );
    if ( x < FPMIN ) return GAMMA + log(x);
    if ( x <= -log(EPS) ) {
        sum = 0.0;
        fact = 1.0;
        for ( k = 1; k <= MAXIT; ++k ) {
            fact *= x/k;
            term = fact/k;
            sum += term;
            if ( term < EPS*sum ) break;
        }
        if ( k > MAXIT ) nrerror( "Series failed in ei." );
        return GAMMA + sum + log(x);
    } else {
        sum = 0.0;
        fact = 1.0;
        for ( k = 1; k <= MAXIT; ++k ) {
            prev = term;
            term *= k/x;
            if ( term < EPS ) break;
            if ( term < prev ) sum += term;
            else {
                sum -= prev;
                break;
            }
        }
        return (1.0 + sum)*exp(x)/x;
    }
}



