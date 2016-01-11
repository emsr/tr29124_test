
#include <stdio.h>
#include <math.h>

#include "nric.h"

int main() {

    dvector3 E1 = { 1.0, 0.0, 0.0 }, E2 = { 0.0, 1.0, 0.0 }, E3 = { 0.0, 0.0, 1.0 };
    dvector3 W = { 0.5, 1.5, 2.5 }, X = { 2.5, -1.5, 0.5 }, Y = { 1.5, -4.5, 1.5 }, Z = { -2.5, 6.5, 5.5 };
    dvector3 X2 = { 5.0, -3.0, 1.0 };
    dmatrix3 A, B, C = { 3.0, 4.5, -1.5, 6.0, 3.5, -0.5, -2.5, -4.5, 5.5 }, D;

    printf( "\n\n  This is a test of the 3 and 2 dimensional vector space structs and thier suppoerting routines.\n" );

    printf( "\n  E1 = %6.3f %6.3f %6.3f", E1.x, E1.y, E1.z );
    printf( "\n  E2 = %6.3f %6.3f %6.3f", E2.x, E2.y, E2.z );
    printf( "\n  E3 = %6.3f %6.3f %6.3f", E3.x, E3.y, E3.z );

    printf( "\n  X = %6.3f %6.3f %6.3f", X.x, X.y, X.z );
    printf( "\n  Y = %6.3f %6.3f %6.3f", Y.x, Y.y, Y.z );
    printf( "\n  Z = %6.3f %6.3f %6.3f", Z.x, Z.y, Z.z );

    printf( "\n\n  E1.X = %6.3f", dot_dv3( E1, X ) );
    printf( "\n\n  E2.X = %6.3f", dot_dv3( E2, X ) );
    printf( "\n\n  E3.X = %6.3f", dot_dv3( E3, X ) );
    printf( "\n\n  Y.Z = %6.3f", dot_dv3( Y, Z ) );

    printf( "\n\n  X^Y = %6.3f %6.3f %6.3f", cross_dv3( X, Y ) );
    printf( "\n\n  X^2X = %6.3f %6.3f %6.3f", cross_dv3( X, X2 ) );

    printf( "\n\n           %6.3f %6.3f %6.3f\n  A = XY = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", A = cart_dv3( X, Y ) );
    printf( "\n\n  detA = %6.3f", det_dv3( A ) );

    printf( "\n\n           %6.3f %6.3f %6.3f\n  B = YZ = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", B = cart_dv3( Y, Z ) );
    printf( "\n\n  detB = %6.3f", det_dv3( B ) );
    printf( "\n\n  trB  = %6.3f", trace_dv3( B ) );
    printf( "\n\n           %6.3f %6.3f %6.3f\n  B-1    = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", inv_dv3( B ) );
    printf( "\n\n           %6.3f %6.3f %6.3f\n  B~     = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", trans_dv3( B ) );

    printf( "\n\n           %6.3f %6.3f %6.3f\n  C      = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", C );
    printf( "\n\n  detC = %6.3f", det_dv3( C ) );
    printf( "\n\n  trC  = %6.3f", trace_dv3( C ) );
    printf( "\n\n           %6.3f %6.3f %6.3f\n  C-1    = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", inv_dv3( C ) );
    printf( "\n\n           %6.3f %6.3f %6.3f\n  C~     = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", trans_dv3( C ) );
    printf( "\n\n           %6.3f %6.3f %6.3f\n  C-1C   = %6.3f %6.3f %6.3f\n           %6.3f %6.3f %6.3f", matmat_dv3( inv_dv3( C ), C ) );

    printf( "\n\n" );
}


