/*  Sample program for demonstrating NR linear algebra functions.
*/


#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/* NR include file - turn parameter checking ON */
#define NR_CHECK
#include "nrlin.h"

/* useful macros */

#define True 1
#define False 0
#define X 1
#define Y 2
#define Z 3
#define U 4
#define PIG (3.141592654)
#define DEGREES(x) ((x)/(PIG)*180) // radians to degrees
#define RAD(x) ((x)/180.0*(PIG)) // degrees to radians


void main( void )
{
	int i, cols, rows, error, errors;
	FILE *inpfile;		/* data input  */
	FILE *outfile;		/* data output */
	
	double **F1, **F1copy, **F2, **F1F2;  /* matrix pointers */
	double *a, *b;							/* vector pointers */

	outfile = fopen("nrdemo.out","w");
	inpfile = fopen("nrdemo.dat","r");
	
	F1 = dmatrix( 1,5,1,5 );
	F1copy = dmatrix( 1,5,1,5 );
	F2 = dmatrix( 1,5,1,2 );
	F1F2 = dmatrix( 1,5, 1,2 );
	a = dvector( 1,2 );
	b = dvector( 1,5 );
	
	errors = 0;
	printf("Reading F1\n");
	dReadMatrix( inpfile, F1, &rows, &cols, &error);
	errors += error;
	if ( rows != 5 || cols != 5 ) errors++;
	
	printf("Reading F2\n");
	dReadMatrix( inpfile, F2, &rows, &cols, &error);
	errors != error;
	if ( rows != 5 || cols != 2 ) errors++;
	
	printf("Reading a\n");
	dReadVector( inpfile, a, &rows, &error);
	errors != error;
	if ( rows != 5 ) errors++;
	
	if ( errors > 0 ) {
		printf("Data errors\n");
		exit(1);
	}
	

	dmcopy(F1, 5, 5, F1copy);
	dinverse(F1copy, 5, F1);	/* F1 is now inverted */
	dmdump( outfile, "F1 inverted", F1, 5, 5, "%11.6lf");
	
	dmmult(F1, 5, 5, F2, 5, 2, F1F2);
	
	dmvmult( F1F2, 5, 2, a, 2, b );
	
	dvdump( outfile, "Result b", b, 5, "%11.4lf");
	dvdump( stdout, "Result b", b, 5, "%11.4lf");  /* on screen */
	
	fprintf( outfile, "Memory used for matrices & vectors: %ld\n",
			mem_used() );
	fprintf( stdout, "Memory used for matrices & vectors: %ld\n",
			mem_used() );
	
	free_dvector( b, 1, 5 );
	free_dvector( a, 1, 2 );
	free_dmatrix( F1, 1,5,1,5 );
	free_dmatrix( F1copy, 1,5,1,5 );
	free_dmatrix( F2, 1,5,1,2 );
	free_dmatrix( F1F2, 1,5, 1,2 );

	fclose(inpfile); 
	fclose(outfile);
}


