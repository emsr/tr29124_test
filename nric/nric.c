
#include "nric.h"
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
    Numerical Recipies standard error handler.
*******************************************************************************/

void  nrerror( const char *error_text ) {

    fprintf( stderr, "\n" );
    fprintf( stderr, "\n  Warning!!!  Lark's Vomit!!!" );
    fprintf( stderr, "\n  Numerical Recipes run-time error..." );
    fprintf( stderr, "\n  %s", error_text );
    fprintf( stderr, "\n  . . . now exiting to system . . ." );
    fprintf( stderr, "\n" );
    exit( 1 );
}



void  nrwarning( const char *warning_text ) {

    fprintf( stderr, "\n" );
    fprintf( stderr, "\n  Warning!!!  Lark's Vomit!!!" );
    fprintf( stderr, "\n  Numerical Recipes run-time warning..." );
    fprintf( stderr, "\n  %s", warning_text );
    fprintf( stderr, "\n  . . . execution continuing . . ." );
    fprintf( stderr, "\n" );
}

