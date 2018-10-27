/* nrcheck.c **************************************************************

  NR Library Package - Wrapper functions for thorough checking

  James Trevelyan,  University of Western Australia
  Revision 2   January 1996

  To use this to add safty checks to your program, add the
  module NRCHECK.C to your project, and #define NR_CHECK before
  including "nrlin.h" in your source file.

**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrlin.h"

/* #define TEST  */  /* use to activate test program at end of source */

#define valid_dvector(v,nl,nh) ((v!=(double *)NULL) && \
					((v[(nl)-1]==-322.0)&&(v[(nh)+1]==-722.0)))
#define valid_vector(v,nl,nh) ((v!=(float *)NULL) && \
					((v[(nl)-1]==-322.0)&&(v[(nh)+1]==-722.0)))
#define valid_dvector_b(v) ( (v!=(double *)NULL) && (*(v)==-322.0) )
#define valid_vector_b(v)  ( (v!=(float *)NULL) && (*(v)==-322.0) )

#define valid_matrix(m,nrl,nrh,ncl,nch) ((m!=(float **)NULL) && \
		  ((m[(nrl)-2]==(float*)0x555) \
		  &&(m[(nrl)][(ncl)-1]==-422.00)&&(m[nrh][nch+1]==-822.0)))
#define valid_dmatrix(m,nrl,nrh,ncl,nch) ((m!=(double **)NULL) && \
		  ((m[(nrl)-2]==(double*)0x555) \
		  &&(m[(nrl)][(ncl)-1]==-422.00)&&(m[nrh][nch+1]==-822.0)))

#define valid_matrix_b(m) ( (m != (float **)NULL ) && ((*(m-1)==(float*)0x555)&&(*m[1]==-422.00)))
#define valid_dmatrix_b(m) ((m != (double **)NULL ) && ((*(m-1)==(double*)0x555)&&(*m[1]==-422.00)))

#define valid_d4vector( y ) valid_dvector( y, 1, 4 )
#define valid_d4vector_b( y ) valid_dvector_b( y )
#define valid_dframe( y )   valid_dmatrix( y, 1, 4, 1, 4 )

#define NoText (char *)NULL
#define NoFile (FILE *)NULL

/* Global storage */

static linenum = 0;
static char *filetext = NoText;
static char *funcname = NoText;
static error_connect = 0;
static char ms[50];

/* function declaration */
void nrerror_handler( void(*handler)(char error_text[]) ); /* nrlin.c */

/* Error handler derived from Numerical Recipes standard error handler */
void nrerror_DNR(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	if ( linenum != 0 ){
		fprintf(stderr,"Call to %s at line %d in %s\n",
								 funcname, linenum, filetext);
	} else {
		fprintf(stderr,"Call location unknown\n");
	}
	fprintf(stderr,"Press <ENTER> to return to system...");
	getchar();
	exit(1);
}

/* Wrapper functions call this to save information about caller */
void RegisterCall_DNR( int lineref, char *fileref, char *func )
{
	if ( !error_connect ) nrerror_handler(nrerror_DNR);
	error_connect = 1;
	linenum = lineref;
	filetext = fileref;
	funcname = func;
}

/* Wrapper functions call this to remove information about caller */
void UnRegisterCall_DNR( void )
{
	linenum = 0;
}

/* Check out a text string parameter */
void Text_DNR( char *text, int argnum )
{
	int i, n;
	unsigned char *textp;

	if ( text == NoText ) {
		sprintf( ms, "NULL text pointer (arg %d)", argnum);
		nrerror_DNR( ms );
	} else {
		n = strlen(text);
		if ( n > 255 ) {
			sprintf( ms, "text string > 255 chars (arg %d)", argnum);
			nrerror_DNR( ms );
		} else {
			if ( n > 10 ) n = 10;
			textp = (unsigned char *)text;
			for (i=0; i<n; i++) {
				if ( textp[i] < ' ' || textp[i] > 127 ) {
					if ( textp[i] != '\n' && textp[i] != '\t' ) {
						sprintf( ms, "text string contains strange chars (arg %d)",
								 argnum);
						nrerror_DNR( ms );
					}
				}
			}
		}
	}
}

/* Check out a 4*4 dmatrix parameter */
void dmat4_DNR( double **a, int argnum )
{
	if ( !valid_dmatrix( a, 1, 4, 1, 4 ) ) {
		sprintf( ms, "arg %d is not a valid 4*4 dmatrix", argnum);
		nrerror_DNR( ms );
	}
}

/* Check out a dmatrix parameter */
void dmat_DNR( double **a, int argnum )
{
	if ( !valid_dmatrix_b( a ) ) {
		sprintf( ms, "arg %d is not a valid dmatrix", argnum);
		nrerror_DNR( ms );
	}
}

/* Check out a dvector parameter */
void dvec_DNR( double *a, int argnum )
{
	if ( !valid_dvector_b( a ) ) {
		sprintf( ms, "arg %d is not a valid dvector", argnum);
		nrerror_DNR( ms );
	}
}

void mat_DNR( float **a, int argnum )
{
	if ( !valid_matrix_b( a ) ) {
		sprintf( ms, "arg %d is not a valid matrix", argnum);
		nrerror_DNR( ms );
	}
}

void vec_DNR( float *a, int argnum )
{
	if ( !valid_vector_b( a ) ) {
		sprintf( ms, "arg %d is not a valid vector", argnum);
		nrerror_DNR( ms );
	}
}

/* Check out a file parameter */
void File_DNR( FILE *fileptr, int argnum )
{
	int i;
	if ( fileptr == NoFile ) {
		sprintf( ms, "NULL file pointer (arg %d)", argnum);
		nrerror_DNR( ms );
	} else {
		i = ferror(fileptr);
		if (i) {
			sprintf( ms, "file error %d(?) (arg %d)", argnum);
			nrerror_DNR( ms );
		}
	}
}



/*  WRAPPER FUNCTIONS BEGIN HERE  */

void free_dmatrix_DNR( int lineref, char *fileref,
					double **a, int nrl, int nrh, int ncl, int nch )
{
	RegisterCall_DNR( lineref, fileref, "free_dmatrix" );
	dmat_DNR( a, 1);
	free_dmatrix( a, nrl, nrh, ncl, nch );
	UnRegisterCall_DNR();
}

void free_dvector_DNR( int lineref, char *fileref,
					double *a, int nrl, int nrh )
{
	RegisterCall_DNR( lineref, fileref, "free_dvector" );
	dvec_DNR( a, 1);
	free_dvector( a, nrl, nrh );
	UnRegisterCall_DNR();
}

double **dmatrix_DNR( int lineref, char *fileref,
					int nrl, int nrh, int ncl, int nch )
{
	double **y;

	RegisterCall_DNR( lineref, fileref, "dmatrix" );
	y = dmatrix( nrl, nrh, ncl, nch );
	UnRegisterCall_DNR();
	return(y);
}

double *dvector_DNR( int lineref, char *fileref,
					int ncl, int nch )
{
	double *y;

	RegisterCall_DNR( lineref, fileref, "dvector" );
	y = dvector( ncl, nch );
	UnRegisterCall_DNR();
	return(y);
}

void dmdump_DNR( int lineref, char *fileref,
					FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format)
{
	RegisterCall_DNR( lineref, fileref, "dmdump" );
	File_DNR( outf, 1 );
	Text_DNR( text, 2 );
	Text_DNR( format, 6);
	dmat_DNR( a, 3);
	dmdump( outf, text, a, a_rows, a_cols, format);
	UnRegisterCall_DNR();
}

void dWriteMatrix_DNR( int lineref, char *fileref,
					FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format)
{
	RegisterCall_DNR( lineref, fileref, "dWriteMatrix" );
	File_DNR( outf, 1 );
	Text_DNR( text, 2 );
	Text_DNR( format, 6);
	dmat_DNR( a, 3);
	dWriteMatrix( outf, text, a, a_rows, a_cols, format);
	UnRegisterCall_DNR();
}

void dvdump_DNR( int lineref, char *fileref,
					FILE *outf, char *text, double *a, int a_els, char *format)
{
	RegisterCall_DNR( lineref, fileref, "dvdump" );
	File_DNR( outf, 1 );
	Text_DNR( text, 2 );
	Text_DNR( format, 5);
	dvec_DNR( a, 3);
	dvdump( outf, text, a, a_els, format);
	UnRegisterCall_DNR();
}

void dWriteVector_DNR( int lineref, char *fileref,
					FILE *outf, char *text, double *a, int a_els, char *format)
{
	RegisterCall_DNR( lineref, fileref, "dWriteVector" );
	File_DNR( outf, 1 );
	Text_DNR( text, 2 );
	Text_DNR( format, 5);
	dvec_DNR( a, 3);
	dWriteVector( outf, text, a, a_els, format);
	UnRegisterCall_DNR();
}

void dWriteScalar_DNR( int lineref, char *fileref,
					FILE *outf, char *text, double a, char *format)
{
	RegisterCall_DNR( lineref, fileref, "dWriteScalar" );
	File_DNR( outf, 1 );
	Text_DNR( text, 2 );
	Text_DNR( format, 4);
	dWriteScalar( outf, text, a, format);
	UnRegisterCall_DNR();
}

void dReadMatrix_DNR( int lineref, char *fileref,
				  FILE *datafile, double **y, int *rows, int *cols, int *error)
/* read a matrix _DNR( int lineref, char *fileref,
				  NR type) */
{
	RegisterCall_DNR( lineref, fileref, "dReadMatrix" );
	File_DNR( datafile, 1 );
	dmat_DNR( y, 2);
	dReadMatrix( datafile, y, rows, cols, error);
	UnRegisterCall_DNR();
}

void dReadVector_DNR( int lineref, char *fileref,
				  FILE *datafile, double *y, int *els, int *error)
/* read a vector (NR type) */
{
	RegisterCall_DNR( lineref, fileref, "dReadVector" );
	File_DNR( datafile, 1 );
	dvec_DNR( y, 2);
	dReadVector( datafile, y, els, error);
	UnRegisterCall_DNR();
}

double dReadScalar_DNR( int lineref, char *fileref,
				  FILE *datafile, int *error)
{
	double y;
	RegisterCall_DNR( lineref, fileref, "dReadScalar" );
	File_DNR( datafile, 1 );
	y = dReadScalar( datafile, error);
	UnRegisterCall_DNR();
	return(y);
}

void dmmult_DNR( int lineref, char *fileref,
				 double **a, int a_rows, int a_cols,
				 double **b, int b_rows, int b_cols, double **y)
/* multiply two matrices a, b, result in y. y must not be same as a or b */
{
	RegisterCall_DNR( lineref, fileref, "dmmult" );
	dmat_DNR( a, 1 );
	dmat_DNR( b, 4 );
	dmat_DNR( y, 7 );
	dmmult( a, a_rows, a_cols, b, b_rows, b_cols, y );
	UnRegisterCall_DNR();
}

void dmvmult_DNR( int lineref, char *fileref,
				  double **a, int a_rows, int a_cols, double *b, int b_els, double *y)
/* multiply a matrix a by vector b, result in y. y can be same as b */
{
	RegisterCall_DNR( lineref, fileref, "dmvmult" );
	dmat_DNR( a, 1 );
	dvec_DNR( b, 4 );
	dvec_DNR( y, 6 );
	dmvmult( a, a_rows, a_cols, b, b_els, y );
	UnRegisterCall_DNR();
}

void dmadd_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols, double **b, double **y)
/* add two matrices a, b, result in y. y can be same as a or b */
{
	RegisterCall_DNR( lineref, fileref, "dmadd" );
	dmat_DNR( a, 1 );
	dmat_DNR( b, 4 );
	dmat_DNR( y, 5 );
	dmadd( a, a_rows, a_cols, b, y );
	UnRegisterCall_DNR();
}

void dmsmy_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols, double r, double **y)
/* multiply a by scalar r, result in y. y can be same as a */
{
	RegisterCall_DNR( lineref, fileref, "dmsmy" );
	dmat_DNR( a, 1 );
	dmat_DNR( y, 5 );
	dmsmy( a, a_rows, a_cols, r, y );
	UnRegisterCall_DNR();
}

void dmsub_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols, double **b, double **y)
/* subtract two matrices a, b, result in y. y can be same as a or b */
{
	RegisterCall_DNR( lineref, fileref, "dmsub" );
	dmat_DNR( a, 1 );
	dmat_DNR( b, 4 );
	dmat_DNR( y, 5 );
	dmsub( a, a_rows, a_cols, b, y );
	UnRegisterCall_DNR();
}

void dmtranspose_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols, double **y)
/* transpose matrix a, result in y. y must not be same as a */
{
	RegisterCall_DNR( lineref, fileref, "dmtranspose" );
	dmat_DNR( a, 1 );
	dmat_DNR( y, 4 );
	dmtranspose( a, a_rows, a_cols, y );
	UnRegisterCall_DNR();
}

void dmfillUT_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols)
/* fill upper half of a (square) (lower triangular form) to make a symmetric
	if not square, will do the best possible */
{
	RegisterCall_DNR( lineref, fileref, "dmfillUT" );
	dmat_DNR( a, 1 );
	dmfillUT( a, a_rows, a_cols);
	UnRegisterCall_DNR();
}

void dvadd_DNR( int lineref, char *fileref,
					double *a, int a_els, double *b, double *y)
{
	RegisterCall_DNR( lineref, fileref, "dvadd" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 3 );
	dvec_DNR( y, 4 );
	dvadd( a, a_els, b, y);
	UnRegisterCall_DNR();
}

void dvsub_DNR( int lineref, char *fileref,
					double *a, int a_els, double *b, double *y)
{
	RegisterCall_DNR( lineref, fileref, "dvsub" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 3 );
	dvec_DNR( y, 4 );
	dvsub( a, a_els, b, y);
	UnRegisterCall_DNR();
}

double dvdot_DNR( int lineref, char *fileref,
					double *a, int a_els, double *b)
{
	double y;

	RegisterCall_DNR( lineref, fileref, "dvdot" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 3 );
	y = dvdot( a, a_els, b);
	UnRegisterCall_DNR();
	return(y);
}

double dvmag_DNR( int lineref, char *fileref,
					double *a, int a_els)
{
	double y;

	RegisterCall_DNR( lineref, fileref, "dvmag" );
	dvec_DNR( a, 1 );
	y = dvmag( a, a_els);
	UnRegisterCall_DNR();
	return(y);
}

double dvmnorm_DNR( int lineref, char *fileref,
					double *a, int a_els)
{
	double y;

	RegisterCall_DNR( lineref, fileref, "dvmnorm" );
	dvec_DNR( a, 1 );
	y = dvmnorm( a, a_els);
	UnRegisterCall_DNR();
	return(y);
}

void dvsmy_DNR( int lineref, char *fileref,
					double *a, int a_els, double r, double *y)
{
	RegisterCall_DNR( lineref, fileref, "dvsmy" );
	dvec_DNR( a, 1 );
	dvec_DNR( y, 4 );
	dvsmy( a, a_els, r, y);
	UnRegisterCall_DNR();
}

void dvpiv_DNR( int lineref, char *fileref,
					double *a, int a_els, double r, double *b, double *y)
{
	RegisterCall_DNR( lineref, fileref, "dvpiv" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 4 );
	dvec_DNR( y, 5 );
	dvpiv( a, a_els, r, b, y);
	UnRegisterCall_DNR();
}



void dmcopy_DNR( int lineref, char *fileref,
					double **a, int a_rows, int a_cols, double **b) /* copy a matrix */
{
	RegisterCall_DNR( lineref, fileref, "dmcopy" );
	dmat_DNR( a, 1 );
	dmat_DNR( b, 4 );
	dmcopy( a, a_rows, a_cols,  b);
	UnRegisterCall_DNR();
}

void dvcopy_DNR( int lineref, char *fileref,
					double *a, int a_els, double *y) /* copy a vector */
{
	RegisterCall_DNR( lineref, fileref, "dvcopy" );
	dvec_DNR( a, 1 );
	dvec_DNR( y, 3 );
	dvcopy( a, a_els, y);
	UnRegisterCall_DNR();
}

double dvcomp_DNR( int lineref, char *fileref,
					double *a, int a_els, double *b) /* compare two vectors */
{
	double y;

	RegisterCall_DNR( lineref, fileref, "dvcomp" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 3 );
	y = dvcomp( a, a_els, b);
	UnRegisterCall_DNR();
	return(y);
}

void dgetcolumn_DNR( int lineref, char *fileref,
				  double **G, int col, double *v, int nels) /* get a column from a matrix */
{
	RegisterCall_DNR( lineref, fileref, "dgetcolumn" );
	dmat_DNR( G, 1 );
	dvec_DNR( v, 3 );
	dgetcolumn( G, col, v, nels );
	UnRegisterCall_DNR();
}

void dputcolumn_DNR( int lineref, char *fileref,
				  double *v, int nels, double **G, int col) /* put a column into a matrix */
{
	RegisterCall_DNR( lineref, fileref, "dputcolumn" );
	dmat_DNR( G, 3 );
	dvec_DNR( v, 1 );
	dputcolumn( v, nels, G, col );
	UnRegisterCall_DNR();
}

void dinverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y )
{
	RegisterCall_DNR( lineref, fileref, "dinverse" );
	dmat_DNR( a, 1 );
	dmat_DNR( y, 3 );
	dinverse( a, n, y );
	UnRegisterCall_DNR();
}

void dinverse_mult_DNR( int lineref, char *fileref,
					double **a, int a_rows, double **b, int b_cols, double **y )
{
	RegisterCall_DNR( lineref, fileref, "dinverse_mult" );
	dmat_DNR( a, 1 );
	dmat_DNR( b, 3 );
	dmat_DNR( y, 5 );
	dinverse_mult( a, a_rows, b, b_cols, y);
	UnRegisterCall_DNR();
}


void dPDSinverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y )
{
	RegisterCall_DNR( lineref, fileref, "dPDSinverse" );
	dmat_DNR( a, 1 );
	dmat_DNR( y, 3 );
	dPDSinverse(a, n, y );
	UnRegisterCall_DNR();
}


void dPDS_L_inverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y )
{
	RegisterCall_DNR( lineref, fileref, "dPDS_L_inverse" );
	dmat_DNR( a, 1 );
	dmat_DNR( y, 3 );
	dPDS_L_inverse( a, n, y );
	UnRegisterCall_DNR();
}


void choldc_DNR( int lineref, char *fileref,
				  float **a, int n, float p[])
{
	RegisterCall_DNR( lineref, fileref, "choldc" );
	mat_DNR( a, 1 );
	vec_DNR( p, 3 );
	choldc( a, n, p );
	UnRegisterCall_DNR();
}


void cholsl_DNR( int lineref, char *fileref,
				  float **a, int n, float p[], float b[], float x[])
{
	RegisterCall_DNR( lineref, fileref, "cholsl" );
	mat_DNR( a, 1 );
	vec_DNR( p, 3 );
	vec_DNR( b, 4 );
	vec_DNR( x, 5 );
	cholsl( a, n, p, b, x );
	UnRegisterCall_DNR();
}


void lubksb_DNR( int lineref, char *fileref,
				  float **a, int n, int *indx, float b[])
{
	RegisterCall_DNR( lineref, fileref, "lubksb" );
	mat_DNR( a, 1 );
	vec_DNR( b, 4 );
	lubksb( a, n, indx, b );
	UnRegisterCall_DNR();
}


void ludcmp_DNR( int lineref, char *fileref,
				  float **a, int n, int *indx, float *d)
{
	RegisterCall_DNR( lineref, fileref, "ludcmp" );
	mat_DNR( a, 1 );
	vec_DNR( d, 4 );
	ludcmp( a, n, indx, d );
	UnRegisterCall_DNR();
}



void dcholdc_DNR( int lineref, char *fileref,
				  double **a, int n, double p[])
{
	RegisterCall_DNR( lineref, fileref, "dcholdc" );
	dmat_DNR( a, 3 );
	dvec_DNR( p, 3 );
	dcholdc( a, n, p );
	UnRegisterCall_DNR();
}


void dcholsl_DNR( int lineref, char *fileref,
				  double **a, int n, double p[], double b[], double x[])
{
	RegisterCall_DNR( lineref, fileref, "dcholsl" );
	dmat_DNR( a, 1 );
	dvec_DNR( p, 3 );
	dvec_DNR( b, 4 );
	dvec_DNR( x, 5 );
	dcholsl( a, n, p, b, x );
	UnRegisterCall_DNR();
}


void dlubksb_DNR( int lineref, char *fileref,
				  double **a, int n, int *indx, double b[])
{
	RegisterCall_DNR( lineref, fileref, "dlubksb" );
	dmat_DNR( a, 1 );
	dvec_DNR( b, 4 );
	dlubksb( a, n, indx, b );
	UnRegisterCall_DNR();
}


void dludcmp_DNR( int lineref, char *fileref,
				  double **a, int n, int *indx, double *d)
{
	RegisterCall_DNR( lineref, fileref, "dludcmp" );
	dmat_DNR( a, 1 );
	dvec_DNR( d, 4 );
	dludcmp( a, n, indx, d );
	UnRegisterCall_DNR();
}



double **dframe_DNR( int lineref, char *fileref	)
{
	double **y;

	RegisterCall_DNR( lineref, fileref, "dframe" );
	y = dframe();
	UnRegisterCall_DNR();
	return y;
}


double *d4vector_DNR( int lineref, char *fileref )
{
	double *y;

	RegisterCall_DNR( lineref, fileref, "d4vector" );
	y = d4vector();
	UnRegisterCall_DNR();
	return y;
}


void free_d4vector_DNR( int lineref, char *fileref, double *y )
{
	RegisterCall_DNR( lineref, fileref, "free_d4vector" );
	free_d4vector(y);
	UnRegisterCall_DNR();
}

void free_dframe_DNR( int lineref, char *fileref, double **y )
{
	RegisterCall_DNR( lineref, fileref, "free_d4vector" );
	free_dframe(y);
	UnRegisterCall_DNR();
}


void chainmult_DNR( int lineref, char *fileref,
					double **a, double **b )
{
	RegisterCall_DNR( lineref, fileref, "chainmult" );
	dmat4_DNR( a, 1 );
	dmat4_DNR( b, 2 );
	chainmult( a, b );
	UnRegisterCall_DNR();
}


void rotframe_DNR( int lineref, char *fileref,
					double **a, int axis, double theta)
{
	RegisterCall_DNR( lineref, fileref, "rotframe" );
	dmat4_DNR( a, 1 );
	rotframe( a, axis, theta );
	UnRegisterCall_DNR();
}


void transframe_DNR( int lineref, char *fileref,
					double **a, double dx, double dy, double dz)
{
	RegisterCall_DNR( lineref, fileref, "transframe" );
	dmat4_DNR( a, 1 );
	transframe( a, dx, dy, dz );
	UnRegisterCall_DNR();
}


void orthog_DNR( int lineref, char *fileref,
					double **a )
{
	RegisterCall_DNR( lineref, fileref, "orthog" );
	dmat4_DNR( a, 1 );
	orthog( a );
	UnRegisterCall_DNR();
}


void orthogKI_DNR( int lineref, char *fileref,
					double **a )
{
	RegisterCall_DNR( lineref, fileref, "orthogKI" );
	dmat4_DNR( a, 1 );
	orthogKI( a );
	UnRegisterCall_DNR();
}


void dframeinverse_DNR( int lineref, char *fileref,
					double **E, double **Ein )
{
	RegisterCall_DNR( lineref, fileref, "dframeinverse" );
	dmat4_DNR( E, 1 );
	dmat4_DNR( Ein, 2 );
	dframeinverse( E, Ein );
	UnRegisterCall_DNR();
}


void dcross4_DNR( int lineref, char *fileref,
					double *a, double *b, double *c )
{
	RegisterCall_DNR( lineref, fileref, "dcross4" );
	dvec_DNR( a, 1 );
	dvec_DNR( b, 2 );
	dvec_DNR( c, 3 );
	dcross4( a, b, c );
	UnRegisterCall_DNR();
}


double column_dot_DNR( int lineref, char *fileref,
					double **G1, int col1, double **G2, int col2 )
{
	double y;

	RegisterCall_DNR( lineref, fileref, "dcolumn_dot" );
	dmat4_DNR( G1, 1 );
	dmat4_DNR( G2, 3 );
	y = column_dot( G1, col1, G2, col2 );
	UnRegisterCall_DNR();
	return (y);
}


void make_d_omega_DNR( int lineref, char *fileref,
					double **G1, double **G2, double *w)
{
	RegisterCall_DNR( lineref, fileref, "make_d_omega" );
	dmat4_DNR( G1, 1 );
	dmat4_DNR( G2, 2 );
	make_d_omega( G1, G2, w );
	UnRegisterCall_DNR();
}


void compare_link_frames_DNR( int lineref, char *fileref,
					double **G1, double **G2, double *wp, double *wr, double *ww )
{
	RegisterCall_DNR( lineref, fileref, "compare_link_frames" );
	dmat4_DNR( G1, 1 );
	dmat4_DNR( G2, 2 );
	compare_link_frames( G1, G2, wp, wr, ww );
	UnRegisterCall_DNR();
}


void compare_frames_DNR( int lineref, char *fileref,
				  double **G1, double **G2, double *wp, double *wr, double *ww )
{
	RegisterCall_DNR( lineref, fileref, "compare_frames" );
	dmat4_DNR( G1, 1 );
	dmat4_DNR( G2, 2 );
	compare_frames( G1, G2, wp, wr, ww );
	UnRegisterCall_DNR();
}


void FindAxis_DNR( int lineref, char *fileref,
				  double **T1, double **T2, double *axis,
					 double *sphi, double *cphi,int *twitching)
{
	RegisterCall_DNR( lineref, fileref, "FindAxis" );
	dmat4_DNR( T1, 1 );
	dmat4_DNR( T2, 2 );
	dvec_DNR( axis, 3 );
	FindAxis( T1, T2, axis, sphi, cphi, twitching );
	UnRegisterCall_DNR();
}


void RotateVector_DNR( int lineref, char *fileref,
				  double *v, double *axis, double sphi, double cphi, double *y)
{
	RegisterCall_DNR( lineref, fileref, "RotateVector" );
	dvec_DNR( v, 1 );
	dvec_DNR( axis, 2 );
	dvec_DNR( y, 5 );
	RotateVector( v, axis, sphi, cphi, y );
	UnRegisterCall_DNR();
}


void RotateFrame_DNR( int lineref, char *fileref,
				  double **T1, double *axis, double sphi, double cphi, double **T2)
{
	RegisterCall_DNR( lineref, fileref, "RotateFrame" );
	dmat_DNR( T1, 1 );
	dvec_DNR( axis, 2 );
	dmat_DNR( T2, 5 );
	RotateFrame( T1, axis, sphi, cphi, T2 );
	UnRegisterCall_DNR();
}


#ifdef TEST
#include "nrcheck.h"
void main( void )
{
	/* test program for above utility routines */

	double **a, **b, **c, **bT;
	double *x, *y, *z;
	FILE *infile, *outfile;
	int a_rows, a_cols, b_rows, b_cols, errors, xn, yn;
	char L1[] = {"Matrix A"};
	char form[] = {"%11.4lf"};

	infile = fopen("mat.dat", "r");
	outfile = fopen("mat.out", "w");

	a = dmatrix(1,10,1,10);
	b = dmatrix(1,10,1,10);
	x = dvector(1,10);
	y = dvector(1,10);

	dReadMatrix( infile, a, &a_rows, &a_cols, &errors);
	printf("Reading a done - hit <return> to go on\n");
	getchar();
	dReadMatrix( infile, b, &b_rows, &b_cols, &errors);
	dReadVector( infile, x, &xn, &errors);
	dReadVector( infile, y, &yn, &errors);
	printf("Reading done - hit <return> to go on\n");
	getchar();

	dmdump( stderr, L1, a, a_rows, a_cols, form);
	dmdump( outfile, "Matrix B", b, b_rows, b_cols, "%8.2lf");
	dvdump( outfile, "Vector x", x, xn, form);
	dvdump( outfile, "Vector y", y, yn, "%8.2lf");
	z = dvector( 1, xn );
	dvadd( x, xn, y, z );
	dvdump( outfile, "x + y", z, xn, "%8.2lf");
	dvsub( x, xn, y, z );
	dvdump( outfile, "x - y", z, xn, "%8.2lf");
	dvsmy( x, xn, 2.0, z );
	dvdump( outfile, "2x", z, xn, "%8.2lf");
	printf("Magnitude of 2x: %7.2lf\n", dvmag( z, xn ));
	printf("dot product x.y: %7.2lf\n", dvdot( x, xn, y));

	dmvmult( a, a_rows, a_cols, x, xn, z );
	dvdump( outfile, "Ax", z, xn, "%8.2lf");
	if ( !valid_dmatrix( a, 1, 10, 1, 10) )nrerror("a invalid");
	c = dmatrix( 1, a_rows, 1, b_cols );
	bT = dmatrix( 1, b_cols, 1, b_rows );
	dmtranspose( b, b_rows, b_cols, bT);
	dmdump( outfile, "Matrix B (transposed)", bT, b_cols, b_rows, "%8.2lf");
	dmmult( a, a_rows, a_cols, bT, b_cols, b_rows, c);
	dmdump( outfile, "Matrix AB", c, a_rows, b_rows, "%8.2lf");
	dmdump( outfile, L1, a, a_rows, a_cols, form);

	/*  dmfillUT( a, a_rows, a_cols );
		 dmdump( outfile, "Symmetrified matrix A", a, a_rows, a_cols, "%8.2lf"); */

	reportmemory(stdout);
	free_dmatrix( a, 1, 10, 1, 10);
	free_dmatrix( b, 1, 10, 1, 10);
	free_dmatrix( c, 1, a_rows, 1, b_cols);
   free_dmatrix( bT, 1, b_cols, 1, b_rows);
	free_dvector( x, 1, 10 );
	free_dvector( y, 1, 10 );
	free_dvector( z, 1, xn );

	reportmemory(stdout);
	printf("Type <return> to exit\n");
	getchar();
}
#endif
