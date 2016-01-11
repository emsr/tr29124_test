
#include <math.h>

#include "nric.h"


/*
 *    Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1;
 *    or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
 *    if isign is input as -1.  data is a complex array of length nn (or equivalently, 
 *    a real array of length 2*nn).  nn must be an integer power of 2!!!.
 */
void fourier1( double *data, unsigned long nn, int isign ) {

    unsigned long n, mmax, m, j, istep, i;
    double tempr, tempi, wtemp, wr, wpr, wpi, wi, theta;

    n = nn << 1;

    /*
     *    This is the bit reversal section of the routine.
     */
    j = 1;
    for ( i = 1; i < n; i += 2 ) {
        if ( j > i ) {
            dswap( &data[j], &data[i] );
            dswap( &data[j+1], &data[i+1]);
        }
        m = n >> 1;
        while ( m >= 2 && j > m ) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    /*
     *    Here begins the Danielson-Lanczos section of the routine.
     */
    mmax = 2;
    while ( n > mmax ) {
        istep = mmax << 1;
        theta = isign*PI2/mmax;
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for ( m = 1; m < mmax; m += 2 ) {
            for ( i = m; i <= n; i += istep ) {
                j = i + mmax;
                tempr = wr*data[j] - wi*data[j+1];
                tempi = wr*data[j+1] + wi*data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp=wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
    }
}


void four_real( double *data, unsigned long n, int isign ) {

    unsigned long i, i1, i2, i3, i4, np3;
    double c1 = 0.5, c2, h1r, h1i, h2r, h2i;
    double wr, wi, wpr, wpi, wtemp, theta;

    theta = PI2/(double) (n>>1);
    if ( isign == 1 ) {
        c2 = -0.5;
        four_real( data, n>>1, 1 );
    } else {
        c2 = 0.5;
        theta = -theta;
    }
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    np3 = n + 3;
    for ( i = 2; i <= (n>>2); ++i ) {
        i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1 )));
        h1r =  c1*(data[i1] + data[i3]);
        h1i =  c1*(data[i2] - data[i4]);
        h2r = -c2*(data[i2] + data[i4]);
        h2i =  c2*(data[i1] - data[i3]);
        data[i1] =  h1r + wr*h2r - wi*h2i;
        data[i2] =  h1i + wr*h2i + wi*h2r;
        data[i3] =  h1r - wr*h2r + wi*h2i;
        data[i4] = -h1i + wr*h2i + wi*h2r;
        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
    }
    if ( isign == 1 ) {
        data[1] = (h1r = data[1]) + data[2];
        data[2] = h1r - data[2];
    } else {
        data[1] = c1*((h1r = data[1]) + data[2]);
        data[2] = c1*(h1r - data[2]);
        fourier1( data, n>>1, -1 );
    }
}


/*
 *    Calculates the sine transform of a set of n real-valued data points stored
 *    in the array y[1..n].  The number n  must be a power of 2.  On exit y is
 *    replaced by its transform.   This routine also calculates the inverse transform
 *    times n/2.
 */
void four_sin( double *y, int n ) {

    int j, np2 = n+2;
    double sum, y1, y2;
    double theta, wi = 0.0, wr = 1.0, wpi, wpr, wtemp;

    theta = PI/n;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    y[1] = 0.0;
    for ( j = 2; j <= (n>>1)+1; ++j ) {
        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
        y1 = wi*(y[j] + y[np2-j]);
        y2 = 0.5*(y[j] - y[np2-j]);
        y[j] = y1 + y2;
        y[np2-j] = y1 - y2;
    }
    four_real( y, n, 1 );
    y[1] *= 0.5;
    sum = y[2] = 0.0;
    for ( j = 1; j <= n-1; j += 2 ) {
        sum += y[j];
        y[j] = y[j+1];
        y[j+1] = sum;
    }
}


/*
 *    Calculates the cosine transform of a set of n real-valued data points stored
 *    in the array y[1..n].  The number n  must be a power of 2.  On exit y is
 *    replaced by its transform.   This routine also calculates the inverse transform
 *    times n/2.
 */
void four_cos1( double *y, int n ) {

    int j, np2 = n+2;
    double sum, y1, y2;
    double theta, wi = 0.0, wr = 1.0, wpi, wpr, wtemp;

    theta = PI/n;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    sum = 0.5*(y[1] - y[n+1]);
    y[1] = 0.5*(y[1] + y[n+1]);
    for ( j = 2; j <= (n>>1); ++j ) {
        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
        y1 = 0.5*(y[j] + y[np2-j]);
        y2 = (y[j] - y[np2-j]);
        y[j] = y1 - wi*y2;
        y[np2-j] = y1  + wi*y2;
        sum += wr*y2;
    }
    four_real( y, n, 1 );
    y[n+1] = y[2];
    y[2] = sum;
    for ( j = 4; j <= n; j += 2 ) {
        sum += y[j];
        y[j] = sum;
    }
}



void four_cos2( double *y, int n ) {


}




