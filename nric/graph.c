
#include "nric.h"
#include <stdio.h>
#include <math.h>

double fresnel_C( double x ) { double c, s; fresnel( x, &c, &s ); return c; }
double fresnel_S( double x ) { double c, s; fresnel( x, &c, &s ); return s; }

double ci( double x ) { double c, s; cisi( x, &c, &s ); return c; }
double si( double x ) { double c, s; cisi( x, &c, &s ); return s; }

/*
double fock_r( double x ) { return fock( x ).r; }
double fock_i( double x ) { return fock( x ).i; }
*/

double E_( double x ) { return legendre_e( x, 1.0 ); }
double F( double x ) { return legendre_f( x, 1.0 ); }
double Pi( double x ) { return legendre_pi( x, 1, 1.0 ); }

double sn( double x ) {
    double sn, cn, dn;
    jacobian_sncndn( x, 0.5, &sn, &cn, &dn );
    return sn;
}

double cn( double x ) {
    double sn, cn, dn;
    jacobian_sncndn( x, 0.5, &sn, &cn, &dn );
    return cn;
}

double dn( double x ) {
    double sn, cn, dn;
    jacobian_sncndn( x, 0.5, &sn, &cn, &dn );
    return dn;
}

double ai( double x ) {
    double ai, bi, aip, bip;
    airy( x, &ai, &bi, &aip, &bip );
    return ai;
}

double bi( double x ) {
    double ai, bi, aip, bip;
    airy( x, &ai, &bi, &aip, &bip );
    return bi;
}

double aip( double x ) {
    double ai, bi, aip, bip;
    airy( x, &ai, &bi, &aip, &bip );
    return aip;
}

double bip( double x ) {
    double ai, bi, aip, bip;
    airy( x, &ai, &bi, &aip, &bip );
    return bip;
}

double bessel_new( double x ) {
    double rj, ry, rjp, ryp;
    bessel_jy( x, 1.0, &rj, &ry, &rjp, &ryp );
    return rj;
}

double chebyshev_t2( double x ) { return chebyshev_t( 2, x ); }
double chebyshev_t3( double x ) { return chebyshev_t( 3, x ); }
double chebyshev_t4( double x ) { return chebyshev_t( 4, x ); }
double chebyshev_t5( double x ) { return chebyshev_t( 5, x ); }

double bessel_j2( double x ) { return bessel_j( 2, x ); }
double bessel_j3( double x ) { return bessel_j( 3, x ); }
double bessel_j4( double x ) { return bessel_j( 4, x ); }
double bessel_j5( double x ) { return bessel_j( 5, x ); }

double bessel_y2( double x ) { return bessel_y( 2, x ); }
double bessel_y3( double x ) { return bessel_y( 3, x ); }
double bessel_y4( double x ) { return bessel_y( 4, x ); }
double bessel_y5( double x ) { return bessel_y( 5, x ); }

double bessel_i2( double x ) { return bessel_i( 2, x ); }
double bessel_i3( double x ) { return bessel_i( 3, x ); }
double bessel_i4( double x ) { return bessel_i( 4, x ); }
double bessel_i5( double x ) { return bessel_i( 5, x ); }

double bessel_k2( double x ) { return bessel_k( 2, x ); }
double bessel_k3( double x ) { return bessel_k( 3, x ); }
double bessel_k4( double x ) { return bessel_k( 4, x ); }
double bessel_k5( double x ) { return bessel_k( 5, x ); }

double hermite_h2( double x ) { return hermite_h( 2, x ); }
double hermite_h3( double x ) { return hermite_h( 3, x ); }
double hermite_h4( double x ) { return hermite_h( 4, x ); }
double hermite_h5( double x ) { return hermite_h( 5, x ); }

double legendre_p2( double x ) { return legendre_p( 2, x ); }
double legendre_p3( double x ) { return legendre_p( 3, x ); }
double legendre_p4( double x ) { return legendre_p( 4, x ); }
double legendre_p5( double x ) { return legendre_p( 5, x ); }

double legendre_q2( double x ) { return legendre_q( 2, x ); }
double legendre_q3( double x ) { return legendre_q( 3, x ); }
double legendre_q4( double x ) { return legendre_q( 4, x ); }
double legendre_q5( double x ) { return legendre_q( 5, x ); }

double sph_bessel_j2( double x ) { return sph_bessel_j( 2, x ); }
double sph_bessel_j3( double x ) { return sph_bessel_j( 3, x ); }
double sph_bessel_j4( double x ) { return sph_bessel_j( 4, x ); }
double sph_bessel_j5( double x ) { return sph_bessel_j( 5, x ); }

double sph_bessel_y2( double x ) { return sph_bessel_y( 2, x ); }
double sph_bessel_y3( double x ) { return sph_bessel_y( 3, x ); }
double sph_bessel_y4( double x ) { return sph_bessel_y( 4, x ); }
double sph_bessel_y5( double x ) { return sph_bessel_y( 5, x ); }

double sph_bessel_i2( double x ) { return sph_bessel_i( 2, x ); }
double sph_bessel_i3( double x ) { return sph_bessel_i( 3, x ); }
double sph_bessel_i4( double x ) { return sph_bessel_i( 4, x ); }
double sph_bessel_i5( double x ) { return sph_bessel_i( 5, x ); }

double sph_bessel_k2( double x ) { return sph_bessel_k( 2, x ); }
double sph_bessel_k3( double x ) { return sph_bessel_k( 3, x ); }
double sph_bessel_k4( double x ) { return sph_bessel_k( 4, x ); }
double sph_bessel_k5( double x ) { return sph_bessel_k( 5, x ); }

double exp_int0( double x ) { return exp_int( 0, x ); }
double exp_int1( double x ) { return exp_int( 1, x ); }
double exp_int2( double x ) { return exp_int( 2, x ); }
double exp_int3( double x ) { return exp_int( 3, x ); }
double exp_int4( double x ) { return exp_int( 4, x ); }
double exp_int5( double x ) { return exp_int( 5, x ); }



int main() {

    double x1;
    double x2;
    char crap[2];
    char title1[40+1];
    char title2[40+1];
    char titlex[40+1];
    char titley[100+1];

/*
    while ( 1 ) {

        printf( "\n Enter x1 x2 (x1 = x2 to stop): " );
        scanf( "%lf %lf", &x1, &x2 );
        if ( x1 == x2 ) break;

        gets( crap );

        printf( "\n Enter title1: " );
        gets( title1 );

        printf( "\n Enter title2: " );
        gets( title2 );

        printf( "\n Enter x title: " );
        gets( titlex );

        printf( "\n Enter y title: " );
        gets( titley );
    }
*/

    plot_func( E_, 0.001, PI/2.0-0.001, "Legendre Elliptic Finction E", "k = 1", "x", "E" );
    plot_func( F, 0.001, PI/2.0-0.001, "Legendre Elliptic Finction F", "k = 1", "x", "F" );
    plot_func( Pi, 0.001, PI/2.0-0.001, "Legendre Elliptic Finction Pi", "k = 1", "x", "Pi" );

    plot_func( sn, -5.0, 5.0, "Jacobian Elliptic Finction sn", "", "x", "sn" );
    plot_func( cn, -5.0, 5.0, "Jacobian Elliptic Finction cn", "", "x", "cn" );
    plot_func( dn, -5.0, 5.0, "Jacobian Elliptic Finction dn", "", "x", "dn" );

    plot_func( ai, -5.0, 5.0, "Airy Finction Ai", "", "x", "Ai" );
    plot_func( bi, -5.0, 5.0, "Airy Finction Bi", "", "x", "Bi" );
    plot_func( aip, -5.0, 5.0, "Airy Finction Ai'", "", "x", "Ai'" );
    plot_func( bip, -5.0, 5.0, "Airy Finction Bi'", "", "x", "Bi'" );
    plot_func( bessel_new, 0.001, 10.0, "Bessel Function J_1", "Order 1", "x", "J1" );


    plot_func( bessel_j0, 0.0, 10.0, "Bessel Function J_0", "Order 0", "x", "J0" );
    plot_func( bessel_j1, 0.0, 10.0, "Bessel Function J_1", "Order 1", "x", "J1" );
    plot_func( bessel_j2, 0.0, 10.0, "Bessel Function J_2", "Order 2", "x", "J2" );
    plot_func( bessel_j3, 0.0, 10.0, "Bessel Function J_3", "Order 3", "x", "J3" );
    plot_func( bessel_j4, 0.0, 10.0, "Bessel Function J_4", "Order 4", "x", "J4" );
    plot_func( bessel_j5, 0.0, 10.0, "Bessel Function J_5", "Order 5", "x", "J5" );

    plot_func( bessel_y0, 0.5, 10.0, "Bessel Function Y_0", "Order 0", "x", "Y0" );
    plot_func( bessel_y1, 0.5, 10.0, "Bessel Function Y_1", "Order 1", "x", "Y1" );
    plot_func( bessel_y2, 0.5, 10.0, "Bessel Function Y_2", "Order 2", "x", "Y2" );
    plot_func( bessel_y3, 0.5, 10.0, "Bessel Function Y_3", "Order 3", "x", "Y3" );
    plot_func( bessel_y4, 0.5, 10.0, "Bessel Function Y_4", "Order 4", "x", "Y4" );
    plot_func( bessel_y5, 0.5, 10.0, "Bessel Function Y_5", "Order 5", "x", "Y5" );

    plot_func( bessel_i0, 0.0, 10.0, "Modified Bessel Function I_0", "Order 0", "x", "I0" );
    plot_func( bessel_i1, 0.0, 10.0, "Modified Bessel Function I_1", "Order 1", "x", "I1" );
    plot_func( bessel_i2, 0.0, 10.0, "Modified Bessel Function I_2", "Order 2", "x", "I2" );
    plot_func( bessel_i3, 0.0, 10.0, "Modified Bessel Function I_3", "Order 3", "x", "I3" );
    plot_func( bessel_i4, 0.0, 10.0, "Modified Bessel Function I_4", "Order 4", "x", "I4" );
    plot_func( bessel_i5, 0.0, 10.0, "Modified Bessel Function I_5", "Order 5", "x", "I5" );

    plot_func( bessel_k0, 0.5, 10.0, "Modified Bessel Function K_0", "Order 0", "x", "K0" );
    plot_func( bessel_k1, 0.5, 10.0, "Modified Bessel Function K_1", "Order 1", "x", "K1" );
    plot_func( bessel_k2, 0.5, 10.0, "Modified Bessel Function K_2", "Order 2", "x", "K2" );
    plot_func( bessel_k3, 0.5, 10.0, "Modified Bessel Function K_3", "Order 3", "x", "K3" );
    plot_func( bessel_k4, 0.5, 10.0, "Modified Bessel Function K_4", "Order 4", "x", "K4" );
    plot_func( bessel_k5, 0.5, 10.0, "Modified Bessel Function K_5", "Order 5", "x", "K5" );

    plot_func( hermite_h0, 0.0, 10.0, "Hermite Polynomial H_0", "Order 0", "x", "H0" );
    plot_func( hermite_h1, 0.0, 10.0, "Hermite Polynomial H_1", "Order 1", "x", "H1" );
    plot_func( hermite_h2, 0.0, 10.0, "Hermite Polynomial H_2", "Order 2", "x", "H2" );
    plot_func( hermite_h3, 0.0, 10.0, "Hermite Polynomial H_3", "Order 3", "x", "H3" );
    plot_func( hermite_h4, 0.0, 10.0, "Hermite Polynomial H_4", "Order 4", "x", "H4" );
    plot_func( hermite_h5, 0.0, 10.0, "Hermite Polynomial H_5", "Order 5", "x", "H5" );

    plot_func( ln_gamma, 0.5, 10.0, "Log Gamma Function", "", "x", "Y0" );

    plot_func( legendre_p0, -1.0, +1.0, "Legendre Polynomial", "Order 0", "x", "P0" );
    plot_func( legendre_p1, -1.0, +1.0, "Legendre Polynomial", "Order 1", "x", "P1" );
    plot_func( legendre_p2, -1.0, +1.0, "Legendre Polynomial", "Order 2", "x", "P2" );
    plot_func( legendre_p3, -1.0, +1.0, "Legendre Polynomial", "Order 3", "x", "P3" );
    plot_func( legendre_p4, -1.0, +1.0, "Legendre Polynomial", "Order 4", "x", "P4" );
    plot_func( legendre_p5, -1.0, +1.0, "Legendre Polynomial", "Order 5", "x", "P5" );

    plot_func( legendre_q0, -0.99, +0.99, "Legendre Polynomial", "Order 0", "x", "Q0" );
    plot_func( legendre_q1, -0.99, +0.99, "Legendre Polynomial", "Order 1", "x", "Q1" );
    plot_func( legendre_q2, -0.99, +0.99, "Legendre Polynomial", "Order 2", "x", "Q2" );
    plot_func( legendre_q3, -0.99, +0.99, "Legendre Polynomial", "Order 3", "x", "Q3" );
    plot_func( legendre_q4, -0.99, +0.99, "Legendre Polynomial", "Order 4", "x", "Q4" );
    plot_func( legendre_q5, -0.99, +0.99, "Legendre Polynomial", "Order 5", "x", "Q5" );

    plot_func( ci, 0.2, 20.0, "Cosine Integral", "", "x", "Ci" );
    plot_func( si, 0.0, 20.0, "Sine Integral", "", "x", "Si" );

    plot_func( dawson, 0.0, 0.3, "Dawson Integral", "", "x", "D" );
    plot_func( dawson, 0.0, 5.0, "Dawson Integral", "", "x", "D" );

    plot_func( fresnel_c, 0.0, 5.0, "Fresnel Cosine Integral", "", "x", "C" );
    plot_func( fresnel_s, 0.0, 5.0, "Fresnel Sine Integral", "", "x", "S" );
    plot_func( fresnel_C, 0.0, 5.0, "New Fresnel Cosine Integral", "", "x", "C" );
    plot_func( fresnel_S, 0.0, 5.0, "New Fresnel Sine Integral", "", "x", "S" );
    plot_func( fock_r, 0.0, 0.5, "Fock Integral", "Real Part", "x", "ReF" );
    plot_func( fock_i, 0.0, 0.5, "Fock Integral", "Imaginary Part", "x", "ImF" );
    plot_func( fock_r, 0.0, 10.0, "Fock Integral", "Real Part", "x", "ReF" );
    plot_func( fock_i, 0.0, 10.0, "Fock Integral", "Imaginary Part", "x", "ImF" );

    plot_func( sph_bessel_j0, 0.0, 10.0, "Spherical Bessel Function j_0", "Order 0", "x", "j0" );
    plot_func( sph_bessel_j1, 0.0, 10.0, "Spherical Bessel Function j_1", "Order 1", "x", "j1" );
    plot_func( sph_bessel_j2, 0.0, 10.0, "Spherical Bessel Function j_2", "Order 2", "x", "j2" );
    plot_func( sph_bessel_j3, 0.0, 10.0, "Spherical Bessel Function j_3", "Order 3", "x", "j3" );
    plot_func( sph_bessel_j4, 0.0, 10.0, "Spherical Bessel Function j_4", "Order 4", "x", "j4" );
    plot_func( sph_bessel_j5, 0.0, 10.0, "Spherical Bessel Function j_5", "Order 5", "x", "j5" );

    plot_func( sph_bessel_y0, 0.5, 10.0, "Spherical Bessel Function y_0", "Order 0", "x", "y0" );
    plot_func( sph_bessel_y1, 0.5, 10.0, "Spherical Bessel Function y_1", "Order 1", "x", "y1" );
    plot_func( sph_bessel_y2, 0.5, 10.0, "Spherical Bessel Function y_2", "Order 2", "x", "y2" );
    plot_func( sph_bessel_y3, 0.5, 10.0, "Spherical Bessel Function y_3", "Order 3", "x", "y3" );
    plot_func( sph_bessel_y4, 0.5, 10.0, "Spherical Bessel Function y_4", "Order 4", "x", "y4" );
    plot_func( sph_bessel_y5, 0.5, 10.0, "Spherical Bessel Function y_5", "Order 5", "x", "y5" );

    plot_func( sph_bessel_i0, 0.0, 10.0, "Modified Spherical Bessel Function i_0", "Order 0", "x", "i0" );
    plot_func( sph_bessel_i1, 0.0, 10.0, "Modified Spherical Bessel Function i_1", "Order 1", "x", "i1" );
    plot_func( sph_bessel_i2, 0.0, 10.0, "Modified Spherical Bessel Function i_2", "Order 2", "x", "i2" );
    plot_func( sph_bessel_i3, 0.0, 10.0, "Modified Spherical Bessel Function i_3", "Order 3", "x", "i3" );
    plot_func( sph_bessel_i4, 0.0, 10.0, "Modified Spherical Bessel Function i_4", "Order 4", "x", "i4" );
    plot_func( sph_bessel_i5, 0.0, 10.0, "Modified Spherical Bessel Function i_5", "Order 5", "x", "i5" );

    plot_func( sph_bessel_k0, 0.5, 10.0, "Modified Spherical Bessel Function k_0", "Order 0", "x", "k0" );
    plot_func( sph_bessel_k1, 0.5, 10.0, "Modified Spherical Bessel Function k_1", "Order 1", "x", "k1" );
    plot_func( sph_bessel_k2, 0.5, 10.0, "Modified Spherical Bessel Function k_2", "Order 2", "x", "k2" );
    plot_func( sph_bessel_k3, 0.5, 10.0, "Modified Spherical Bessel Function k_3", "Order 3", "x", "k3" );
    plot_func( sph_bessel_k4, 0.5, 10.0, "Modified Spherical Bessel Function k_4", "Order 4", "x", "k4" );
    plot_func( sph_bessel_k5, 0.5, 10.0, "Modified Spherical Bessel Function k_5", "Order 5", "x", "k5" );

    plot_func( chebyshev_t0, -1.0, 1.0, "Chebyshev Polynomial T_0", "Order 0", "x", "T0" );
    plot_func( chebyshev_t1, -1.0, 1.0, "Chebyshev Polynomial T_1", "Order 1", "x", "T1" );
    plot_func( chebyshev_t2, -1.0, 1.0, "Chebyshev Polynomial T_2", "Order 2", "x", "T2" );
    plot_func( chebyshev_t3, -1.0, 1.0, "Chebyshev Polynomial T_3", "Order 3", "x", "T3" );
    plot_func( chebyshev_t4, -1.0, 1.0, "Chebyshev Polynomial T_4", "Order 4", "x", "T4" );
    plot_func( chebyshev_t5, -1.0, 1.0, "Chebyshev Polynomial T_5", "Order 5", "x", "T5" );

    plot_func( exp_int0, 0.1, 5.0, "Exponential Integral Function E_0", "Order 0", "x", "E0" );
    plot_func( exp_int1, 0.1, 5.0, "Exponential Integral Function E_1", "Order 1", "x", "E1" );
    plot_func( exp_int2, 0.1, 5.0, "Exponential Integral Function E_2", "Order 2", "x", "E2" );
    plot_func( exp_int3, 0.1, 5.0, "Exponential Integral Function E_3", "Order 3", "x", "E3" );
    plot_func( exp_int4, 0.1, 5.0, "Exponential Integral Function E_4", "Order 4", "x", "E4" );
    plot_func( exp_int5, 0.1, 5.0, "Exponential Integral Function E_5", "Order 5", "x", "E5" );

    plot_func( ei, 0.1, 5.0, "Exponential Integral Function Ei", "", "x", "Ei" );

}

