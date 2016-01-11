
#include <stdio.h>
#include <math.h>


#include "nric.h"

main() {

    double **a, **ainv, **ain, **b, **bin;
    double *temp31, **temp33, **temp35;
    double **r, **v;
    double **IU, **IV, *w;
    double *d;
    double p, det, trace;
    int *i;
    int j, k, s;

    IU = dmatrix( 1, 3, 1, 3 );
    IV = dmatrix( 1, 3, 1, 3 );
    r = dmatrix( 1, 3, 1, 3 );

    v = dmatrix( 1, 3, 1, 3 );
    w = dvector( 1, 3 );
    i = ivector( 1, 3 );

    ain = dmatrix( 1, 3, 1, 3 );
    a = dmatrix( 1, 3, 1, 3 );
    ainv = dmatrix( 1, 3, 1, 3 );

    bin = dmatrix( 1, 3, 1, 5 );
    b = dmatrix( 1, 3, 1, 5 );

    d = dvector( 1, 3 );

    temp31 = dvector( 1, 3 );
    temp33 = dmatrix( 1, 3, 1, 3 );
    temp35 = dmatrix( 1, 3, 1, 5 );

    ain[1][1] = 1.339673;
    ain[1][2] = 2.240185;
    ain[1][3] = 3.708725;
    ain[2][1] = 2.205473;
    ain[2][2] = 4.991265;
    ain[2][3] = 3.093471;
    ain[3][1] = 7.653201;
    ain[3][2] = 6.592175;
    ain[3][3] = 9.102383;

    bin[1][1] = 3.419463;
    bin[2][1] = 2.204656;
    bin[3][1] = 3.285281;
    bin[1][2] = 4.090252;
    bin[2][2] = 9.937392;
    bin[3][2] = 0.729572;
    bin[1][3] = 2.374950;
    bin[2][3] = 6.374950;
    bin[3][3] = 5.104705;
    bin[1][4] = 8.590235;
    bin[2][4] = 2.495969;
    bin[3][4] = 5.308507;
    bin[1][5] = 7.882022;
    bin[2][5] = 7.703502;
    bin[3][5] = 0.320573;

    dcopymat( ain, 1, 3, 1, 3, a );
    dcopymat( bin, 1, 3, 1, 5, b );


    printf( "\n" );
    printf( "\n Input matrix" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );

    printf( "\n" );
    printf( "\n Input vectors of Gauss-Jordan elimination" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n B = | %9.6f %9.6f %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n",
      b[1][1], b[1][2], b[1][3], b[1][4], b[1][5],
      b[2][1], b[2][2], b[2][3], b[2][4], b[2][5],
      b[3][1], b[3][2], b[3][3], b[3][4], b[3][5] );

    gauss_jordan( a, 3, b, 5 );

    printf( "\n" );
    printf( "\n Output inverse matrix of Gauss-Jordan elimination" );
    printf( "\n" );
    printf( "        | %9.6f %9.6f %9.6f |\n A^-1 = | %9.6f %9.6f %9.6f |\n        | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );

    printf( "\n" );
    printf( "\n Output solution vectors of Gauss-Jordan elimination" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n B = | %9.6f %9.6f %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n",
      b[1][1], b[1][2], b[1][3], b[1][4], b[1][5],
      b[2][1], b[2][2], b[2][3], b[2][4], b[2][5],
      b[3][1], b[3][2], b[3][3], b[3][4], b[3][5] );

    dmatmat( ain, a, 1, 3, 1, 3, 1, 3, temp33 );

    printf( "\n" );
    printf( "\n Verify A.A^-1 = I" );
    printf( "\n" );
    printf( "        | %9.6f %9.6f %9.6f |\n A^-1 = | %9.6f %9.6f %9.6f |\n        | %9.6f %9.6f %9.6f |\n",
      temp33[1][1], temp33[1][2], temp33[1][3],
      temp33[2][1], temp33[2][2], temp33[2][3],
      temp33[3][1], temp33[3][2], temp33[3][3] );

    dmatmat( ain, b, 1, 3, 1, 3, 1, 5, temp35 );

    printf( "\n" );
    printf( "\n Verify A.X = B" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n B = | %9.6f %9.6f %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f %9.6f %9.6f |\n",
      temp35[1][1], temp35[1][2], temp35[1][3], temp35[1][4], temp35[1][5],
      temp35[2][1], temp35[2][2], temp35[2][3], temp35[2][4], temp35[2][5],
      temp35[3][1], temp35[3][2], temp35[3][3], temp35[3][4], temp35[3][5] );



    dcopymat( ain, 1, 3, 1, 3, a );

    printf( "\n" );
    printf( "\n Input matrix" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );


    lu_decomp( a, 3, i, &p );


    printf( "\n" );
    printf( "\n Output matrix of LU decompostion" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );

    printf( "\n" );
    printf( "\n Output row permutation vector of LU decomposition" );
    printf( "\n" );
    printf( "     | %d |\n I = | %d |\n     | %d |\n",
      i[1], i[2], i[3] );

    printf( "\n" );
    printf( "\n Output parity of LU decomposition" );
    printf( "\n" );
    printf( "      P = %9.6f\n", p );

    printf( "\n" );
    printf( "\n Reconstruction of input matrix from LU decomposition" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1],         a[1][2],                         a[1][3],
      a[2][1]*a[1][1], a[2][1]*a[1][2]+a[2][2],         a[2][1]*a[1][3]+a[2][3],
      a[3][1]*a[1][1], a[3][1]*a[1][2]+a[3][2]*a[2][2], a[3][1]*a[1][3]+a[3][2]*a[2][3]+a[3][3] );

    det = lu_determinant( a, 3, p );

    printf( "\n" );
    printf( "\n Determinant of input matrix" );
    printf( "\n" );
    printf( "\n detA = %9.6f\n", det );

    trace = lu_trace( a, 3 );

    printf( "\n" );
    printf( "\n Trace of input matrix" );
    printf( "\n" );
    printf( "\n trA = %9.6f\n", trace );

    lu_invert( a, 3, i, ainv );

    printf( "\n" );
    printf( "\n Inverse of input matrix" );
    printf( "\n" );
    printf( "        | %9.6f %9.6f %9.6f |\n A^-1 = | %9.6f %9.6f %9.6f |\n        | %9.6f %9.6f %9.6f |\n",
      ainv[1][1], ainv[1][2], ainv[1][3],
      ainv[2][1], ainv[2][2], ainv[2][3],
      ainv[3][1], ainv[3][2], ainv[3][3] );

    dmatmat( ainv, ain, 1, 3, 1, 3, 1, 3, temp33 );

    printf( "\n" );
    printf( "\n Verify A^-1.A = I" );
    printf( "\n" );
    printf( "          | %9.6f %9.6f %9.6f |\n A^-1.A = | %9.6f %9.6f %9.6f |\n          | %9.6f %9.6f %9.6f |\n",
      temp33[1][1], temp33[1][2], temp33[1][3],
      temp33[2][1], temp33[2][2], temp33[2][3],
      temp33[3][1], temp33[3][2], temp33[3][3] );



    dcopymat( ain, 1, 3, 1, 3, a );


    printf( "\n" );
    printf( "\n Input matrix" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );


    sv_decomp( a, 3, 3, w, v );


    printf( "\n" );
    printf( "\n Output matrix of SV decompostion" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n U = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );

    printf( "\n" );
    printf( "\n Output vector of SV decompostion" );
    printf( "\n" );
    printf( "     | %9.6f |\n W = | %9.6f |\n     | %9.6f |\n",
      w[1], w[2], w[3] );

    printf( "\n" );
    printf( "\n Output matrix of SV decompostion" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n V = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      v[1][1], v[1][2], v[1][3],
      v[2][1], v[2][2], v[2][3],
      v[3][1], v[3][2], v[3][3] );


    for ( j = 1; j <= 3; ++j ) {
        for ( k = 1; k <= 3; ++k ) {
            r[j][k] = 0.0;
            IU[j][k] = 0.0;
            IV[j][k] = 0.0;
            for ( s = 1; s <= 3; ++s ) {
                r[j][k] += a[j][s]*w[s]*v[k][s];
                IU[j][k] += (a[s][j])*(a[s][k]);
                IV[j][k] += (v[s][j])*(v[s][k]);
            }
        }
    }

    printf( "\n" );
    printf( "\n Verify U~.U = I" );
    printf( "\n" );
    printf( "        | %9.6f %9.6f %9.6f |\n U~.U = | %9.6f %9.6f %9.6f |\n        | %9.6f %9.6f %9.6f |\n",
      IU[1][1], IU[1][2], IU[1][3],
      IU[2][1], IU[2][2], IU[2][3],
      IU[3][1], IU[3][2], IU[3][3] );

    printf( "\n" );
    printf( "\n Verify V~.V = I" );
    printf( "\n" );
    printf( "        | %9.6f %9.6f %9.6f |\n V~.V = | %9.6f %9.6f %9.6f |\n        | %9.6f %9.6f %9.6f |\n",
      IV[1][1], IV[1][2], IV[1][3],
      IV[2][1], IV[2][2], IV[2][3],
      IV[3][1], IV[3][2], IV[3][3] );

    printf( "\n" );
    printf( "\n Reconstruction of input matrix from SV decomposition" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      r[1][1], r[1][2], r[1][3],
      r[2][1], r[2][2], r[2][3],
      r[3][1], r[3][2], r[3][3] );



    dcopymat( ain, 1, 3, 1, 3, a );
    for ( j = 1; j <= 3; ++j ) for ( k = 1; k < j; ++k ) a[j][k] = a[k][j];


    printf( "\n" );
    printf( "\n Input matrix" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );


    cholesky_decomp( a, 3, d );


    printf( "\n" );
    printf( "\n Output matrix of Cholesky decompostion" );
    printf( "\n" );
    printf( "     | %9.6f %9.6f %9.6f |\n A = | %9.6f %9.6f %9.6f |\n     | %9.6f %9.6f %9.6f |\n",
      a[1][1], a[1][2], a[1][3],
      a[2][1], a[2][2], a[2][3],
      a[3][1], a[3][2], a[3][3] );

    printf( "\n" );
    printf( "\n Output vector of Cholesky decompostion" );
    printf( "\n" );
    printf( "     | %9.6f |\n p = | %9.6f |\n     | %9.6f |\n",
      d[1], d[2], d[3] );




    free_dmatrix( IU, 1, 3, 1, 3 );
    free_dmatrix( IV, 1, 3, 1, 3 );
    free_dmatrix( r, 1, 3, 1, 3 );
    free_dmatrix( bin, 1, 3, 1, 5 );
    free_dmatrix( ain, 1, 3, 1, 3 );
    free_dmatrix( temp33, 1, 3, 1, 3 );
    free_dmatrix( temp35, 1, 3, 1, 5 );
    free_dmatrix( b, 1, 3, 1, 5 );
    free_dmatrix( a, 1, 3, 1, 3 );
    free_dmatrix( v, 1, 3, 1, 3 );
    free_dvector( temp31, 1, 3 );
    free_dvector( w, 1, 3 );
    free_dvector( d, 1, 3 );
    free_ivector( i, 1, 3 );
}

