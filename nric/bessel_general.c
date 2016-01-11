


#include <math.h>


#include "nric.h"


void bessel_jy( double xnu, double x, double *rj, double *ry, double *rjp, double *ryp ) {

    const double  EPS    =  1.0e-16;
    const double  FPMIN  =  1.0e-30;
    const int     MAXIT  =  10000;
    const double  XMIN   =  2.0;

    int i, isign, l, nl;
    double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
           fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
           rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
           temp, w, x2, xi, xi2, xmu, xmu2;

    if ( x <= 0.0 || xnu < 0.0 ) nrerror( "Bad arguments in bessel_jy." );

    nl = ( x < XMIN ? (int)(xnu + 0.5) : imax( 0, (int)(xnu - x + 1.5) ) );

    xmu = xnu - nl;
    xmu2 = xmu*xmu;
    xi = 1.0/x;
    xi2 = 2.0*xi;
    w = xi2/PI;
    isign = 1;
    h = xnu*xi;
    if ( h < FPMIN ) h = FPMIN;
    b = xi2*xnu;
    d = 0.0;
    c = h;
    for ( i = 1; i <= MAXIT; ++i ) {
        b += xi2;
        d = b - d;
        if ( fabs(d) < FPMIN ) d = FPMIN;
        c = b - 1.0/c;
        if ( fabs(c) < FPMIN ) c = FPMIN;
        d = 1.0/d;
        del = c*d;
        h = del*h;
        if ( d < 0.0 ) isign = -isign;
        if ( fabs(del - 1.0) < EPS ) break;
    }
    if ( i > MAXIT ) nrerror( "Argument x too large in bessel_jy; try asymptotic expansion." );
    rjl = isign*FPMIN;
    rjpl = h*rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = xnu*xi;
    for ( l = nl; l >= 1; --l ) {
        rjtemp = fact*rjl + rjpl;
        fact -= xi;
        rjpl = fact*rjtemp - rjl;
        rjl = rjtemp;
    }
    if ( rjl == 0.0 ) rjl = EPS;
    f= rjpl/rjl;
    if ( x < XMIN ) {
        x2 = 0.5*x;
        pimu = PI*xmu;
        fact = ( fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu) );
        d = -log(x2);
        e = xmu*d;
        fact2 = ( fabs(e) < EPS ? 1.0 : sinh(e)/e );
        bessel_cheb( xmu, &gam1, &gam2, &gampl, &gammi );
        ff = (2.0/PI)*fact*(gam1*cosh(e) + gam2*fact2*d);
        e = exp(e);
        p = e/(PI*gampl);
        q = 1.0/(e*PI*gammi);
        pimu2 = 0.5*pimu;
        fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2 );
        r = PI*pimu2*fact3*fact3;
        c = 1.0;
        d = -x2*x2;
        sum = ff + r*q;
        sum1 = p;
        for ( i = 1; i <= MAXIT; ++i ) {
            ff = (i*ff + p + q)/(i*i - xmu2);
            c *= d/i;
            p /= i - xmu;
            q /= i + xmu;
            del = c*(ff + r*q);
            sum += del; 
            del1 = c*p - i*del;
            sum1 += del1;
            if ( fabs(del) < EPS*(1.0 + fabs(sum)) ) break;
        }
        if ( i > MAXIT ) nrerror("Bessel y series failed to converge in bessel_jy." );
        rymu = -sum;
        ry1 = -sum1*xi2;
        rymup = xmu*xi*rymu - ry1;
        rjmu = w/(rymup - f*rymu);
    } else {
        a = 0.25 - xmu2;
        q = 1.0;
        p = -0.5*xi;
        br = 2.0*x;
        bi = 2.0;
        fact = a*xi/(p*p + q*q);
        cr = br + q*fact;
        ci = bi + p*fact;
        den = br*br + bi*bi;
        dr = br/den;
        di = -bi/den;
        dlr = cr*dr - ci*di;
        dli = cr*di + ci*dr;
        temp = p*dlr - q*dli;
        q = p*dli + q*dlr;
        p = temp;
        for ( i = 2; i <= MAXIT; ++i ) {
            a += 2*(i - 1);
            bi += 2.0;
            dr = a*dr + br;
            di = a*di + bi;
            if ( fabs(dr) + fabs(di) < FPMIN ) dr = FPMIN;
            fact = a/(cr*cr + ci*ci);
            cr = br + cr*fact;
            ci = bi - ci*fact;
            if ( fabs(cr) + fabs(ci) < FPMIN ) cr = FPMIN;
            den = dr*dr + di*di;
            dr /= den;
            di /= -den;
            dlr = cr*dr - ci*di;
            dli = cr*di + ci*dr;
            temp = p*dlr - q*dli;
            q = p*dli + q*dlr;
            p = temp;
            if ( fabs(dlr - 1.0) + fabs(dli) < EPS ) break;
        }
        if ( i > MAXIT ) nrerror( "Lentz's method failed in bessel_jy." );
        gam = (p - f)/q;
        rjmu = sqrt(w/((p - f)*gam + q));
        rjmu = dsign( rjmu, rjl );
        rymu = rjmu*gam;
        rymup = rymu*(p + q/gam);
        ry1 = xmu*xi*rymu - rymup;
    }
    fact = rjmu/rjl;
    *rj = rjl1*fact;
    *rjp = rjp1*fact;
    for ( i = 1; i <= nl; ++i ) {
        rytemp = (xmu + i)*xi2*ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
    }
    *ry = rymu;
    *ryp = xnu*xi*rymu - ry1;
}



void bessel_ik( double xnu, double x, double *ri, double *rk, double *rip, double *rkp ) {

    const double EPS = 1.0e-16;
    const double FPMIN = 1.0e-30;
    const int MAXIT = 10000;
    const double XMIN = 2.0;

    int i, l, nl;
    double a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, ff,
      gam1, gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, ril1, rimu,
      rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xi, xi2, xmu, xmu2;

    if ( x <= 0.0 || xnu < 0.0 ) nrerror( "Bad arguments in bessel_ik." );

    nl = (int)(xnu + 0.5);

    xmu = xnu - nl;
    xmu2 = xmu*xmu;
    xi = 1.0/x;
    xi2 = 2.0*xi;
    h = xnu*xi;
    if ( h < FPMIN ) h = FPMIN;
    b = xi2*xnu;
    d = 0.0;
    c = h;
    for ( i = 1; i <= MAXIT; ++i ) {
        b += xi2;
        d = 1.0/(b + d);
        c = b + 1.0/c;
        del = c*d;
        h = del*h;
        if ( fabs(del-1.0) < EPS ) break;
    }
    if ( i > MAXIT ) nrerror( "Argument x too large in bessel_ik; try asymptotic expansion." );
    ril = FPMIN;
    ripl = h*ril;
    ril1 = ril;
    rip1 = ripl;
    fact = xnu*xi;
    for ( l = nl; l >= 1; --l ) {
        ritemp = fact*ril + ripl;
        fact -= xi;
        ripl = fact*ritemp + ril;
        ril = ritemp;
    }
    f= ripl/ril;
    if ( x < XMIN ) {
        x2 = 0.5*x;
        pimu = PI*xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        e = xmu*d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
        bessel_cheb( xmu, &gam1, &gam2, &gampl, &gammi );
        ff = fact*(gam1*cosh(e) + gam2*fact2*d);
        sum = ff;
        e = exp(e);
        p = 0.5*e/gampl;
        q = 0.5/(e*gammi);
        c = 1.0;
        d = x2*x2;
        sum1 = p;
        for ( i = 1; i <= MAXIT; ++i ) {
            ff = (i*ff + p + q)/(i*i - xmu2);
            c *= d/i;
            p /= i - xmu;
            q /= i + xmu;
            del = c*ff;
            sum += del; 
            del1 = c*(p - i*ff);
            sum1 += del1;
            if ( fabs(del) < EPS*fabs(sum) ) break;
        }
        if ( i > MAXIT ) nrerror("Bessel k series failed to converge in bessel_jy." );
        rkmu = sum;
        rk1 = sum1*xi2;
    } else {
        b = 2.0*(1.0 + x);
        d = 1.0/b;
        h = delh = d;
        q1 = 0.0;
        q2 = 1.0;
        a1 = 0.25 - xmu2;
        q = c = a1;
        a = -a1;
        s = 1.0 + q*delh;
        for ( i = 2; i <= MAXIT; ++i ) {
            a -= 2*(i - 1);
            c = -a*c/i;
            qnew = (q1 - b*q2)/a;
            q1 = q2;
            q2 = qnew;
            q += c*qnew;
            b += 2.0;
            d = 1.0/(b + a*d);
            delh = (b*d - 1.0)*delh;
            h += delh;
            dels = q*delh;
            s += dels;
            if ( fabs(dels/s) < EPS ) break;
        }
        if ( i > MAXIT ) nrerror( "Steed's method failed in bessel_ik." );
        h = a1*h;
        rkmu = sqrt(PI/(2.0*x))*exp(-x)/s;
        rk1 = rkmu*(xmu + x + 0.5 - h)*xi;
    }
    rkmup = xmu*xi*rkmu - rk1;
    rimu = xi/(f*rkmu - rkmup);
    *ri = rimu*ril1/ril;
    *rip = rimu*rip1/ril;
    for ( i = 1; i <= nl; ++i ) {
        rktemp = (xmu + i)*xi2*rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    *rk = rkmu;
    *rkp = xnu*xi*rkmu - rk1;
}



void bessel_cheb( double x, double *gam1, double *gam2, double *gampl, double *gammi ) {

    double xx;
    static double c1[] = {
      -1.142022680371168e0, 6.5165112670737e-3,
      3.087090173086e-4, -3.4706269649e-6,
      6.9437664e-9, 3.67765e-11, -1.356e-13 };
    static double c2[] = {
      1.843740587300905e0, -7.86528408447867e-2,
      1.2719271366546e-3, -4.9717367042e-6,
      -3.31261198e-8, 2.423096e-10, -1.702e-13, -1.49e-15 };

    xx = 8.0*x*x - 1.0;
    *gam1 = chebyshev_eval( -1.0, +1.0, c1, 7, xx );
    *gam2 = chebyshev_eval( -1.0, +1.0, c2, 8, xx );
    *gampl = *gam2 - x*(*gam1);
    *gammi = *gam2 + x*(*gam1);
}



void airy( double x, double *ai, double *bi, double *aip, double *bip ) {

    double absx, ri, rip, rj, rjp, rk, rkp, rootx, ry, ryp, z;

    absx = fabs(x);
    rootx = sqrt(absx);
    z = 2.0*absx*rootx/3.0;
    if ( x > 0.0 ) {

        bessel_ik( 1.0/3.0, z, &ri, &rk, &rip, &rkp );
        *ai = rootx*rk/(SQRT3*PI);
        *bi = rootx*(rk/PI + 2.0*ri/SQRT3);

        bessel_ik( 2.0/3.0, z, &ri, &rk, &rip, &rkp );
        *aip = -x*rk/(SQRT3*PI);
        *bip = x*(rk/PI + 2.0*ri/SQRT3);

    } else if ( x < 0.0 ) {

        bessel_jy( 1.0/3.0, z, &rj, &ry, &rjp, &ryp );
        *ai = 0.5*rootx*(rj - ry/SQRT3);
        *bi = -0.5*rootx*(ry + rj/SQRT3);

        bessel_jy( 2.0/3.0, z, &rj, &ry, &rjp, &ryp );
        *aip = 0.5*absx*(ry/SQRT3 + rj);
        *bip = 0.5*absx*(rj/SQRT3 - ry);

    } else {

       *ai = 0.35512420385917071;
       *bi = *ai*SQRT3;

       *aip = -0.25874932837133380;
       *bip = -*aip*SQRT3;

    }
}



void wairy( double x, dcomplex *w1, dcomplex *w2, dcomplex *w1p, dcomplex *w2p ) {
    double factor = SQRTPI;
    double ai, bi, aip, bip;

    airy( x, &ai, &bi, &aip, &bip );
    *w1 = mk_dc( ai, bi );
    *w1 = rmul_dc( factor, *w1 );

    *w2 = mk_dc( ai, -bi );
    *w2 = rmul_dc( factor, *w2 );

    *w1p = mk_dc( aip, bip );
    *w1p = rmul_dc( factor, *w1p );

    *w2p = mk_dc( aip, -bip );
    *w2p = rmul_dc( factor, *w2p );

}


void airy_root( double *alpha, double *beta, double *alpha_p, double *beta_p, int s ) {
    const double _alpha[10] = {
         -2.33810741,
         -4.08794944,
         -5.52055983,
         -6.78670809,
         -7.94413359,
         -9.02265085,
        -10.04017434,
        -11.00852430,
        -11.93601556,
        -12.82877675
    };
    const double _beta[10] = {
         -1.17371322,
         -3.27109330,
         -4.83073784,
         -6.16985213,
         -7.37676208,
         -8.49194885,
         -9.53819438,
        -10.52991351,
        -11.47695355,
        -12.38641714
    };
    const double _alpha_p[10] = {
         -1.01879297,
         -3.24819758,
         -4.82009921,
         -6.16330736,
         -7.37217726,
         -8.48848673,
         -9.53544905,
        -10.52766040,
        -11.47505663,
        -12.38478837
    };
    const double _beta_p[10] = {
         -2.29443968,
         -4.07315509,
         -5.51239573,
         -6.78129445,
         -7.94017869,
         -9.01958336,
        -10.03769633,
        -11.00646267,
        -11.93426165,
        -12.82725831
    };

    if ( s < 1 ) {
        nrerror( "invalid root number in airy_root." );
    }
    if ( s > 10 ) {
        return;
    }
    else {
        *alpha = _alpha[s-1];
        *beta = _beta[s-1];
        *alpha_p = _alpha_p[s-1];
        *beta_p = _beta_p[s-1];
    }
}



void sph_bessel( int n, double x, double *sj, double *sy, double *sjp, double *syp ) {

    double factor, order, rj, rjp, ry, ryp;

    if ( n < 0 || x < 0.0 ) nrerror( "Bad arguments in sph_bessel." );

    order = n + 0.5;

    bessel_jy( x, order, &rj, &ry, &rjp, &ryp );
    factor = SQRTPIO2/sqrt(x);
    *sj = factor*rj;
    *sy = factor*ry;
    *sjp = factor*rjp - *sj/(2.0*x);
    *syp = factor*ryp - *sy/(2.0*x);
}


