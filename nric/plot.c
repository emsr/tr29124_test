
#include "nric.h"
#include <stdio.h>
#include <math.h>
#include <string.h>



void plot_func(
  double (*funk)(double),
  double x1,
  double x2,
  const char *t1,
  const char *t2,
  const char *tx,
  const char *ty
)
{


    int jz, j, i;
    double ysml, ybig, x, dyj, dx, *y;
    char **screen;
    char format_labelx[40];
    char *title1, *title2, *titlex, *titley;
    int len1, len2, lenx, leny;

    const int ISCREEN = 100;
    const int JSCREEN = 41;
    const char PLUS = '+';
    const char BLANK = ' ';
    const char ZERO = '-';
    const char YY = '|';
    const char XX = '-';
    const char FF = 'x';

    screen = cmatrix( 1, ISCREEN, 1, JSCREEN );
    y = dvector( 1, ISCREEN );
    title1 = cvector( 1, ISCREEN );
    title2 = cvector( 1, ISCREEN );
    titlex = cvector( 1, ISCREEN );
    titley = cvector( 1, JSCREEN );

    /*
     *    Massage titles.
     */
    len1 = ( strlen( t1 ) < ISCREEN ? strlen( t1 ) : ISCREEN );
    for ( j = 1; j <= ISCREEN; ++j ) { title1[j] = BLANK; }
    for ( j = 1; j <= len1; ++j ) { title1[(int)( ( ISCREEN - len1)/2.0) + j] = t1[j-1]; }

    len2 = ( strlen( t2 ) < ISCREEN ? strlen( t2 ) : ISCREEN );
    for ( j = 1; j <= ISCREEN; ++j ) { title2[j] = BLANK; }
    for ( j = 1; j <= len2; ++j ) { title2[(int)( ( ISCREEN - len2)/2.0) + j] = t2[j-1]; }

    lenx = ( strlen( tx ) < ISCREEN ? strlen( tx ) : ISCREEN );
    for ( j = 1; j <= ISCREEN; ++j ) { titlex[j] = BLANK; }
    for ( j = 1; j <= lenx; ++j ) { titlex[(int)( ( ISCREEN - lenx)/2.0) + j] = tx[j-1]; }

    leny = ( strlen( ty ) < JSCREEN ? strlen( ty ) : JSCREEN );
    for ( j = 1; j <= JSCREEN; ++j ) { titley[j] = BLANK; }
    for ( j = 1; j <= leny; ++j ) { titley[(int)( ( JSCREEN + leny)/2.0) - j] = ty[j-1]; }


    screen[1][1] = screen[1][JSCREEN] = screen[ISCREEN][1] = screen[ISCREEN][JSCREEN] = PLUS;

    for( j = 2; j <= (JSCREEN-1); ++j ) {

        screen[1][j] = screen[ISCREEN][j] = YY;
    }
    for( i = 2; i <= (ISCREEN-1); ++i ) {

        screen[i][1] = screen[i][JSCREEN] = XX;
        for ( j = 2; j <= (JSCREEN-1); ++j ) {

            screen[i][j] = BLANK;
        }
    }
    dx = ( x2 - x1 )/(ISCREEN-1);
    x = x1;

    /*
     *    Evaluate the function over the requested interval.
     *    Keep track of the maximum and minimum values of the function.
     *    With this logic, the interval from ysml to ybig will always
     *    include the x axis (y == 0).
     */
    ysml = ybig = 0.0;
    for ( i = 1; i <= ISCREEN; ++i ) {

        y[i] = (*funk)(x);
        if ( y[i] < ysml ) ysml = y[i];
        if ( y[i] > ybig ) ybig = y[i];
        x += dx;
    }
    if ( ybig == ysml ) ybig = ysml + 1.0;
    dyj = (JSCREEN-1)/( ybig - ysml);
    jz = 1 - (int)( ysml*dyj);
    for ( i = 1; i <= ISCREEN; ++i ) {
        screen[i][jz] = ZERO;
        j = 1 + (int)( ( y[i] - ysml)*dyj);
        screen[i][j] = FF;
    }

    /*
     *    Print title(s)
     */
    if ( len1 ) {
        printf( "\n" );
        if ( leny ) { printf( "%3s", "" ); }
        for ( i = 1; i <= 12; ++i ) { printf( "%c", BLANK ); }
        for ( i = 1; i <= ISCREEN; ++i ) { printf( "%c", title1[i] ); }
        printf( "\n" );
    }
    if ( len2 ) {
        printf( "\n" );
        if ( leny ) { printf( "%3s", "" ); }
        for ( i = 1; i <= 12; ++i ) { printf( "%c", BLANK ); }
        for ( i = 1; i <= ISCREEN; ++i ) { printf( "%c", title2[i] ); }
        printf( "\n" );
    }
    printf( "\n" );

    /*
     *    Print upper limit and top line.
     */
    if ( leny ) { printf( " %c ", titley[JSCREEN] ); }
    printf( " %10.3f ", ybig );
    for ( i = 1; i <= ISCREEN; ++i ) printf( "%c", screen[i][JSCREEN] );
    printf( "\n" );

    /*
     *    Print graph.
     */
    for ( j = (JSCREEN-1); j >= 2; --j ) {
        if ( leny ) { printf( " %c ", titley[j] ); }
        for ( i = 1; i <= 12; ++i ) { printf( "%c", BLANK ); }
        for ( i = 1; i <= ISCREEN; ++i ) printf( "%c", screen[i][j] );
        printf( "\n" );
    }

    /*
     *    Print lower limit and bottom line.
     */
    if ( leny ) { printf( " %c ", titley[1] ); }
    printf( " %10.3f ", ysml );
    for ( i = 1; i <= ISCREEN; ++i ) printf( "%c", screen[i][1] );
    printf( "\n" );

    /*
     *    Print lower and upper x limits.
     */
    if ( leny ) { printf( "%3s", "" ); }
    sprintf( format_labelx, "%%5s %%10.3f %%%ds %%10.3f\n", ISCREEN - 13 );
    printf( format_labelx, "", x1, "", x2 );

    if ( lenx ) {
        printf( "\n" );
        if ( leny ) { printf( "%3s", "" ); }
        for ( i = 1; i <= 12; ++i ) { printf( "%c", BLANK ); }
        for ( i = 1; i <= ISCREEN; ++i ) { printf( "%c", titlex[i] ); }
        printf( "\n" );
    }
    printf( "\n" );
}

