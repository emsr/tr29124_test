
#include <math.h>

#include "nric.h"


/*
 *  Fit data to a straight line. NRiC pp 523-528.
 */

void  fit( double *x, double *y, int ndata, double *sigma, int mwt, 
           double *a, double *b, double *sigmaa,
           double *sigmab, double *chisqr, double *q)  {

    int  i;
    double  wt, t, temp, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat;

    *b = 0.0;
    if ( mwt ) {

        /*
         *    Accumulate sums.
         */
        ss = 0.0;
        for ( i = 1; i <= ndata; i++) {
            wt = 1.0/( sigma[i]*sigma[i] );
            ss += wt;
            sx += x[i]*wt;
            sy += y[i]*wt;
        }
    } else {
        for ( i = 1; i <= ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata;
    }
    sxoss = sx/ss;

    if ( mwt ) {
        for ( i = 1; i <= ndata; i++ ) {
            t = ( x[i] - sxoss)/sigma[i];
            st2 += t*t;
            *b += t*y[i]/sigma[i];
        }
    } else {
        for ( i = 1; i <= ndata; i++ ) {
            t = x[i] - sxoss;
            st2 += t*t;
            *b += t*y[i];
        }
    }
    *b /= st2;
    *a = (sy - sx*(*b))/ss;
    *sigmaa = sqrt((1.0 + sx*sx/(ss*st2))/ss);
    *sigmab = sqrt(1.0/st2);
    *chisqr = 0.0;

    if ( !mwt )  {
        for ( i = 1; i <= ndata; i++ ) *chisqr += (temp = y[i] - *a - *b*x[i], temp*temp);
        *q = 1.0;
        sigdat = sqrt(*chisqr/(ndata - 2));
        *sigmaa *= sigdat;
        *sigmab *= sigdat;
    }
    else  {
        for ( i = 1; i <= ndata; i++ ) *chisqr += (temp = y[i] - *a - *b*x[i], temp*temp)/sigma[i];
        *q = gamma_q( 0.5*(ndata - 2), 0.5*(*chisqr));
    }
}


/*
 *    Given a set of data points x[1..ndata], y[1..ndata] with individual standard deviations
 *    sig[1..ndata], use chi^2 minimization to determine the coefficients a[1..ma] of the
 *    fitting function y = S a_i.afunk_i(x).  The fitting coefficients are solved using
 *    Singular value decomposition of the ndata by ma matrix.  The Arrays u[1..ndata][1..ma]
 *    v[1..ma][1..ma], and w[1..ma] provide workspace on input, on output they define the SVD,
 *    and can be used to determine the covariance matrix.  The input array tol[1..ma] is used as
 *    the tolerance (relative to the largest singular value element wmax) of the singular value
 *    array w[1..ma].  Elements of w[1..ma] such that w[k] < tol[k]*wmax will be set to zero.
 *    The program returns the fit parameters a[1..ma] and the chi squared statistic, chisqr.
 *    The user supplies a routine funcs that returns the ma basis functions evaluated at x
 *    in the array afunc[1..ma].
 */
void svd_fit( double *x, double *y, double *sig, int ndata, double *a, int ma,
              double **u, double **v, double *w, double *tol, double *chisqr,
              void (*funcs)( double, double *, int) ) {

    int i, j;
    double wmax, temp, sum, *b, *afunc;

    b = dvector( 1, ndata );
    afunc = dvector( 1, ma );

    for ( i = 1; i <= ndata; ++i ) {
        (*funcs)( x[i], afunc, ma );
        for ( j = 1; j <= ma; ++j ) u[i][j] = afunc[j]/sig[i];
        b[i] = y[i]/sig[i];
    }

    sv_decomp( u, ndata, ma, w, v );

    /*
     *    Find max singular value and edit singular value vector.
     */
    wmax = 0.0;
    for ( j = 1; j <= ma; ++j ) if ( w[j] > wmax ) wmax = w[j];
    for ( j = 1; j <= ma; ++j ) if ( w[j] < tol[j]*wmax ) w[j] = 0.0;

    sv_backsub( u, w, v, ndata, ma, b, a );

    *chisqr = 0.0;
    for ( i = 1; i <= ndata; ++i ) {
        (*funcs)( x[i], afunc, ma );
        for ( sum = 0.0, j = 1; j <= ma; ++j ) sum += a[j]*afunc[j];
        *chisqr += ( temp = (y[i] - sum)/sig[i], temp*temp );
    }

    free_dvector( afunc, 1, ma );
    free_dvector( b, 1, ndata );
}


/*
 *    Evaluates the covariance matrix cov[1..ma][1..ma] of the fit for ma parameters
 *    obtained by svd_fit.  The input arrays v[1..ma][1..ma] and w[1..ma] are from svd_fit.
 */
void svd_var( double **v, int ma, double *w, double **cov ) {

    int i, j, k;
    double sum, *wt;

    wt = dvector( 1, ma );

    for ( i = 1; i <= ma; ++i ) {
        wt[i] = 0.0;
        if ( w[i] ) wt[i] = 1.0/(w[i]*w[i]);
    }
    for ( i = 1; i <= ma; ++i ) {
        for ( j = 1; j <= i; ++j ) {
            for ( sum = 0.0, k = 1; k <= ma; ++k ) sum += v[i][k]*v[j][k]*wt[k];
            cov[j][i] = cov[i][j] = sum;
        }
    }

    free_dvector( wt, 1, ma );
}
 


void funcs_poly( double x, double *p, int n ) {

    int j;

    p[1] = 1.0;
    for ( j = 2; j <= n; ++j ) p[j] = x*p[j-1];
}


void funcs_legendre( double x, double *p, int n ) {

    int j;
    double twox, f1, f2, d;

    p[1] = 1.0;
    p[2] = x;
    if ( n > 2 ) {
        twox = 2.0*x;
        f2 = x;
        d = 1.0;
        for ( j = 3; j <= n; ++j ) {
            f1 = ++d;
            f2 += twox;
            p[j] = (f2*p[j-1] - f1*p[j-2])/d;
        }
    }
}


/*
 *    Given a set of data points x[1..ndata], y[1..ndata] with individual standard deviations
 *    sig[1..ndata], and an a priori guess for the fitting parameters a[0..maxit][1..ma], use
 *    differential correction to determine the parameters a[0..maxit][1..ma] of the
 *    fitting function y = y(x;a).  Corrections to the fitting parameters are solved using
 *    Singular value decomposition of the ndata by ma matrix.  The Arrays u[1..ndata][1..ma]
 *    v[1..ma][1..ma], and w[1..ma] provide workspace on input, on output they define the SVD,
 *    and can be used to determine the covariance matrix.  The input array tol[1..ma] is used as
 *    the tolerance (relative to the largest singular value element wmax) of the singular value
 *    array w[1..ma].  Elements of w[1..ma] such that w[k] < tol[k]*wmax will be set to zero.
 *    The program returns the fit parameters a[0..maxit][1..ma] and the chi squared statistic,
 *    chisqr[0..maxit].
 *    The user supplies a routine  partials  that returns the computed fitting functions
 *    and the ma partial derivatives of the fitting function with respect to the ma fitting
 *    parameters evaluated at x in the array parts[0..ma].
 */
void svd_diff_corr( double *x, double *y, double *sig, int ndata, double **a, int ma,
                    double **u, double **v, double *w, double *tol, double *chisqr,
                    double **omc, double ***covar, int maxit,
                    void (*partials)( double, double *, int, double * ) ) {

    int it, i, j;
    double wmax, *parts, *da;

    parts = dvector( 1, ma );
    da = dvector( 1, ma );

    for ( it = 0; it <= maxit; ++it ) {

        chisqr[it] = 0.0;
        for ( i = 1; i <= ndata; ++i ) {
            (*partials)( x[i], a[it], ma, parts );
            chisqr[it] += ( omc[it][i] = (y[i] - parts[0])/sig[i], omc[it][i]*omc[it][i] );
            for ( j = 1; j <= ma; ++j ) u[i][j] = parts[j]/sig[i];
        }

        sv_decomp( u, ndata, ma, w, v );

        /*
         *    Find max singular value and edit singular value vector.
         */
        wmax = 0.0;
        for ( j = 1; j <= ma; ++j ) if ( w[j] > wmax ) wmax = w[j];
        for ( j = 1; j <= ma; ++j ) if ( w[j] < tol[j]*wmax ) w[j] = 0.0;

        /*
         *    Compute the covariance.
         */
        svd_var( v, ma, w, covar[it] );

        /*
         *    Backsubstitute for the state vector correction and apply the correction.
         */
        sv_backsub( u, w, v, ndata, ma, omc[it], da );

        if ( it < maxit ) for ( j = 1; j <= ma; ++j ) a[it+1][j] += da[j];

    }

    free_dvector( da, 1, ma );
    free_dvector( parts, 1, ma );
}



