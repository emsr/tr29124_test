

#include <math.h>

#include "nric.h"


double carlson_rf( double x, double y, double z ) {

    double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;
    const double ERRTOL = 0.0025, TINY = 1.5e-38, BIG = 3.0e37,
                 C1 = 1.0/24.0, C2 = 0.1, C3 = 3.0/44.0, C4 = 1.0/14.0;

    if ( dmin3( x, y, z ) < 0.0 || dmin3( x + y, y + z, z + x ) < TINY || dmax3( x, y, z ) > BIG ) nrerror( "Invalid arguments in rf." );

    xt = x;
    yt = y;
    zt = z;
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx*sqrty + sqrty*sqrtz + sqrtz*sqrtx;
        xt = 0.25*(xt + alamb);
        yt = 0.25*(yt + alamb);
        zt = 0.25*(zt + alamb);
        ave = (xt + yt + zt)/3.0;
        delx = (ave - xt)/ave;
        dely = (ave - yt)/ave;
        delz = (ave - zt)/ave;
    } while ( dmax3( fabs(delx), fabs(dely), fabs(delz) ) > ERRTOL );
    e2 = delx*dely - delz*delz;
    e3 = delx*dely*delz;
    return (1.0 + (C1*e2 - C2 - C3*e3)*e2 + C4*e3)/sqrt(ave);
}


double carlson_rd( double x, double y, double z ) {

    double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;
    const double ERRTOL = 0.0015, TINY = 1.0e-25, BIG = 4.5e21,
                 C1 = 3.0/14.0, C2 = 1.0/6.0, C3 = 9.0/22.0, C4 = 3.0/26.0, C5 = 0.25*C3, C6 = 1.5*C4;

    if ( dmin( x, y ) < 0.0 || dmin( x+y, z ) < TINY || dmax3( x, y, z ) > BIG ) nrerror( "Invalid arguments in rd." );

    xt = x;
    yt = y; 
    zt = z;
    sum = 0.0;
    fac = 1.0;
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx*sqrty + sqrty*sqrtz + sqrtz*sqrtx;
        sum += fac/(sqrtz*(zt + alamb));
        fac *= 0.25;
        xt = 0.25*(xt + alamb);
        yt = 0.25*(yt + alamb);
        zt = 0.25*(zt + alamb);
        ave = 0.2*(xt + yt + 3.0*zt);
        delx = (ave - xt)/ave;
        dely = (ave - yt)/ave;
        delz = (ave - zt)/ave;
    } while ( dmax3( fabs(delx), fabs(dely), fabs(delz) ) > ERRTOL );
    ea = delx*dely;
    eb = delz*delz;
    ec = ea - eb;
    ed = ea - 6.0*eb;
    ee = ed + 2.0*ec;
    return 3.0*sum + fac*(1.0 + ed*(-C1 + C5*ed - C6*delz*ee) + delz*(C2*ee + delz*(-C3*ec + delz*C4*ea)))/(ave*sqrt(ave));
}


double carlson_rj( double x, double y, double z, double p ) {

    double a, alamb, alpha, ans, ave, b, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac,
           pt, rcx, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt, zt;
    const double ERRTOL = 0.0015, TINY = 2.5e-13, BIG = 9.0e11, C1 = 3.0/11.0, C2 = 1.0/3.0, C3 = 3.0/22.0, C4 = 3.0/26.0,
                 C5 = 0.75*C3, C6 = 1.5*C4, C7 = 0.5*C2, C8 = 2.0*C3;

    if ( dmin3( x, y, z ) < 0.0 || dmin( dmin3( x+y, y+z, x+x ), fabs(p) ) < TINY || dmax( dmax3( x, y, z ), fabs(p) ) > BIG )
        nrerror( "Invalid arguments in carlson_rj." );

    sum = 0.0;
    fac = 1.0;
    if ( p > 0.0 ) {
        xt = x;
        yt = y;
        zt = z;
    pt = p;
    } else {
        xt = dmin3( x, y, z );
        zt = dmax3( x, y, z );
        yt = x + y + z - xt - zt;
        a = 1.0/(yt - p);
        b = a*(zt - yt)*(yt - xt);
        pt = yt + b;
        rho = xt*zt/yt;
        tau = p*pt/yt;
        rcx = carlson_rc( rho, tau );
    }
    do {
        sqrtx = sqrt(xt);
        sqrty = sqrt(yt);
        sqrtz = sqrt(zt);
        alamb = sqrtx*sqrty + sqrty*sqrtz + sqrtz*sqrtx;
        alpha = dsqr( pt*(sqrtx + sqrty + sqrtz) + sqrtx*sqrty*sqrtz );
        beta = pt*dsqr( pt + alamb );
        sum += fac*carlson_rc( alpha, beta );
        fac *= 0.25;
        xt = 0.25*(xt + alamb);
        yt = 0.25*(yt + alamb);
        zt = 0.25*(zt + alamb);
        pt = 0.25*(pt + alamb);
        ave = 0.2*(xt + yt + zt + pt + pt);
        delx = (ave - xt)/ave;
        dely = (ave - yt)/ave;
        delz = (ave - zt)/ave;
        delp = (ave - pt)/ave;
    } while ( dmax( dmax3( fabs(delx), fabs(dely), fabs(delz) ), fabs(delp) ) > ERRTOL );
    ea = delx*dely + dely*delz + delz*delx;
    eb = delx*dely*delz;
    ec = delp*delp;
    ed = ea - 3.0*ec;
    ee = eb + 2.0*delp*(ea - ec);
    ans = 3.0*sum + fac*(1.0 + ed*(-C1 + C5*ed - C6*ee) + eb*(C7 + delp*(-C8 + delp*C4)) + delp*ea*(C2 - delp*C3) - C2*delp*ec)/(ave*sqrt(ave));
    if ( p < 0.0 ) ans = a*(b*ans + 3.0*(rcx - carlson_rf( xt, yt, zt )));
    return ans;
}


double carlson_rc( double x, double y ) {

    double alamb, ave, s, w, xt, yt;
    const double ERRTOL = 0.0012, TINY = 1.69e-38, SQRTNY = 1.3e-19, BIG = 3.0e37, TNBG = TINY*BIG,
                 COMP1 = 2.236/SQRTNY, COMP2 = TNBG*TNBG/25.0, C1 = 0.3, C2 = 1.0/7.0, C3 = 0.375, C4 = 9.0/22.0;

    if ( x < 0.0 || y == 0.0 || (x + fabs(y)) < TINY || (x + fabs(y)) > BIG || (y < -COMP1 && x > 0.0 && x < COMP2) )
        nrerror( "Invalid arguments in carlson_rc" );

    if ( y > 0.0 ) {
        xt = x;
        yt = y;
        w = 1.0;
    } else {
        xt = x - y;
        yt = -y;
        w = sqrt(x)/sqrt(xt);
    }
    do {
        alamb = 2.0*sqrt(xt)*sqrt(yt) + yt;
        xt = 0.25*(xt + alamb);
        yt = 0.25*(yt + alamb);
        ave = (xt + yt + yt)/3.0;
        s = (yt - ave)/ave;
    } while ( fabs(s) > ERRTOL );
    return w*(1.0 + s*s*(C1 + s*(C2 + s*(C3 + s*C4))))/sqrt(ave);
}


/*
 *    Legendre elliptic function of the second kind, E(phi, k) evaluated using
 *    Carlson's elliptic functions Rf and Rd.  The argument ranges are  0 <= phi <= pi/2
 *    and  0 <= ksin(phi) <= 1.
 */
double legendre_e( double phi, double k ) {

    double s, c, cc, q, ks;

    if ( phi < 0.0 || phi > PI/2.0 ) nrerror( "Argument  phi  out of range in legendre_f." );

    s = sin(phi);
    ks = k*s;
    if ( ks < 0.0 || ks > 1.0 ) nrerror( "Argument  k sin(phi)  out of range in legendre_f." );
    c = cos(phi);
    cc = c*c;
    q = (1.0 - ks)*(1.0 + ks);

    return s*(carlson_rf( cc, q, 1.0 ) - dsqr(ks))*carlson_rd( cc, q, 1.0 );
}


/*
 *    Legendre elliptic function of the first kind, F(phi, k) evaluated using
 *    Carlson's elliptic functions Rf.  The argument ranges are  0 <= phi <= pi/2
 *    and  0 <= ksin(phi) <= 1.
 */
double legendre_f( double phi, double k ) {

    double s, c, cc, q, ks;

    if ( phi < 0.0 || phi > PI/2.0 ) nrerror( "Argument  phi  out of range in legendre_f." );

    s = sin(phi);
    ks = k*s;
    if ( ks < 0.0 || ks > 1.0 ) nrerror( "Argument  k sin(phi)  out of range in legendre_f." );
    c = cos(phi);
    cc = c*c;
    q = (1.0 - ks)*(1.0 + ks);

    return s*carlson_rf( cc, q, 1.0 );
}


/*
 *    Legendre elliptic function of the third kind, Pi(phi,n,k) evaluated using
 *    Carlson's elliptic functions Rj and Rf.  The sign convention on n is opposite 
 *    Abramowtz and Stegun.  The ranges are  0 <= phi <= pi/2  and  0 <= ksin(phi) <= 1.
 */
double legendre_pi( double phi, double n, double k ) {

    double s, c, cc, nss, q, ks;

    if ( phi < 0.0 || phi > PI/2.0 ) nrerror( "Argument  phi  out of range in legendre_f." );

    s = sin(phi);
    ks = k*s;
    if ( ks < 0.0 || ks > 1.0 ) nrerror( "Argument  k sin(phi)  out of range in legendre_f." );
    nss = n*s*s;
    c = cos(phi);
    cc = c*c;
    q = (1.0 - ks)*(1.0 + ks);
    return s*(carlson_rf( cc, q, 1.0 ) - nss*carlson_rj( cc, q, 1.0, 1.0 + nss )/3.0);
}


void jacobian_sncndn( double uu, double mmc, double *sn, double *cn, double *dn ) {

    double a, b, c, d, mc, u;
    double *m, *n;
    int i, ii, l, bo;

    const int N = 20;
    const double CA = 0.00001;

    m = dvector( 1, N );
    n = dvector( 1, N );

    mc = mmc;
    u = uu;
    if ( mc ) {
        bo = (mc < 0.0);
        if ( bo ) {
            d = 1.0 - mc;
            mc /= -1.0/d;
            u *= (d = sqrt(d));
        }
        a = 1.0;
        *dn = 1.0;
        for ( i = 1; i <= N; ++i ) {
            l = i;
            m[i] = a;
            n[i] = (mc = sqrt(mc));
            c = 0.5*(a + mc);
            if ( fabs(a - mc) <= CA*a ) break;
            mc *= a;
            a = c;
        }
        u *= c;
        *sn = sin(u);
        *cn = cos(u);
        if ( *sn ) {
            a = (*cn)/(*sn);
            c *= a;
            for ( ii = l; ii >= 1; --ii ) {
                b = m[ii];
                a *= c;
                c *= (*dn);
                *dn = (n[ii] + a)/(b + a);
                a = c/b;
            }
            a = 1.0/sqrt(c*c + 1.0);
            *sn = dsign( a, *sn );
            *cn = c*(*sn);
        }
        if ( bo ) {
            a = *dn;
            *dn = *cn;
            *cn = a;
            *sn /= d;
        }
    } else {
        *cn = 1.0/cosh(u);
        *dn = *cn;
        *sn = tanh(u);
    }

    free_dvector( n, 1, N );
    free_dvector( m, 1, N );
}


