/* nrutil.c ***************************************************************

  NR Library Package - Utilities Module

  James Trevelyan,  University of Western Australia
  Revision 2   January 1996

**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrlin.h"

/* define V_CHECK for validity checks on matrices to be included */
#define V_CHECK
/* define this if you want to use the test program below */
#undef TEST


void dmmult( double **a, int a_rows, int a_cols, 
             double **b, int b_rows, int b_cols, double **y)
/* multiply two matrices a, b, result in y. y must not be same as a or b */
{
   int i, j, k;
   double sum;

   if ( a_cols != b_rows ) {
      fprintf(stderr,"a_cols <> b_rows (%d,%d): dmmult\n", a_cols, b_rows);
      exit(1);
   }

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmmult\n");    
   if ( !valid_dmatrix_b( b ) ) 
      nrerror("Invalid 2nd matrix: dmmult\n");    
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid result matrix: dmmult\n"); 
#endif

   /*  getchar();
       dmdump( stdout, "Matrix a", a, a_rows, a_cols, "%8.2lf");
       dmdump( stdout, "Matrix b", b, b_rows, b_cols, "%8.2lf");
       getchar();
   */
   for ( i=1; i<=a_rows; i++ ) 
      for ( j=1; j<=b_cols; j++ ) {
         sum = 0.0;
         for ( k=1; k<=a_cols; k++ ) sum += a[i][k]*b[k][j];
         y[i][j] = sum;
      }
}

void dmvmult( double **a, int a_rows, int a_cols, double *b, int b_els, double *y)
/* multiply a matrix a by vector b, result in y. y can be same as b */
{
   int i, k;
   double sum;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid matrix: dmvmult\n");   
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to vector: dmvmult");
#endif

   for ( i=1; i<=a_rows; i++ ) {
      sum = 0.0;
      for ( k=1; k<=a_cols; k++ ) sum += a[i][k]*b[k];
      y[i] = sum;
   }
}

void dmadd( double **a, int a_rows, int a_cols, double **b, double **y)
/* add two matrices a, b, result in y. y can be same as a or b */
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmadd\n"); 
   if ( !valid_dmatrix_b( b ) ) 
      nrerror("Invalid 2nd matrix: dmadd\n"); 
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid result matrix: dmadd\n");  
#endif

   for ( i=1; i<=a_rows; i++ ) 
      for ( j=1; j<=a_cols; j++ ) {
         y[i][j] = a[i][j] + b[i][j];
      }
}

void dmsmy( double **a, int a_rows, int a_cols, double r, double **y)
/* multiply a by scalar r, result in y. y can be same as a */
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmsmy\n"); 
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid result matrix: dmsmy\n");  
#endif

   for ( i=1; i<=a_rows; i++ ) 
      for ( j=1; j<=a_cols; j++ ) {
         y[i][j] = a[i][j] * r;
      }
}

void dmsub( double **a, int a_rows, int a_cols, double **b, double **y)
/* subtract two matrices a, b, result in y. y can be same as a or b */
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmsub\n"); 
   if ( !valid_dmatrix_b( b ) ) 
      nrerror("Invalid 2nd matrix: dmsub\n"); 
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid result matrix: dmsub\n");  
#endif

   for ( i=1; i<=a_rows; i++ ) 
      for ( j=1; j<=a_cols; j++ ) {
         y[i][j] = a[i][j] - b[i][j];
      }
}

void dmtranspose( double **a, int a_rows, int a_cols, double **y)
/* transpose matrix a, result in y. y must not be same as a */
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmtranspose\n");   
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid result matrix: dmtranspose\n");    
#endif

   for ( i=1; i<=a_rows; i++ ) 
      for ( j=1; j<=a_cols; j++ ) {
         y[j][i] = a[i][j];
      }
}

void dmfillUT( double **a, int a_rows, int a_cols)
/* fill upper half of a (square) (lower triangular form) to make a symmetric 
   if not square, will do the best possible */
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid 1st matrix: dmfillUT\n");  
#endif

   for ( i=1; i<=a_rows; i++ ) 
      for ( j=i; j<=a_cols; j++ ) {
         a[j][i] = a[i][j];
      }
}




void dmdump( FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format)
{
   int i, j, k;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid matrix: dmdump\n");    
#endif

   fprintf( outf, "%s", text);
   for ( i=1; i<=a_rows; i++ ) {
      fprintf( outf, "\n  [%d][1]:", i);
      k=0;
      for ( j=1; j<=a_cols; j++ ) {
         fprintf( outf, format, a[i][j]);
         k = (k+1)%6;
         if ( (k == 0) && (j<a_cols) )fprintf(outf, "\n          ");
      }
   }
   fprintf( outf, "\n");
}   

void dWriteMatrix( FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format)
{
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) ) 
      nrerror("Invalid matrix: dmdump\n");    
#endif
   fprintf( outf, "* %s", text);
   fprintf( outf, "\nM %d %d", a_rows, a_cols );
   for ( i=1; i<=a_rows; i++ ) {
      fprintf( outf, "\n");
      for ( j=1; j<=a_cols; j++ ) {
         fprintf( outf, format, a[i][j]);
      }
   }
   fprintf( outf, "\n");
}

void dvdump( FILE *outf, char *text, double *a, int a_els, char *format)
{
   int j, k;

   if ( !valid_dvector_b( a ) ) {
      fprintf( stderr, "%s: Invalid pointer to vector\n", text);
      return;
   }
   fprintf( outf, "%s\n", text);
   k = 0;
   for ( j=1; j<=a_els; j++ ) {
      fprintf( outf, format, a[j]);
      k = (k+1)%6;
      if ( k == 0 )fprintf( outf, "\n");
   }
   fprintf( outf, "\n");
}   

void dWriteVector( FILE *outf, char *text, double *a, int a_els, char *format)
{
   int j;

   if ( !valid_dvector_b( a ) ) {
      fprintf( stderr, "%s: Invalid pointer to vector\n", text);
      return;
   }
   fprintf( outf, "* %s\n", text);
   fprintf( outf, "V %d\n", a_els);
   for ( j=1; j<=a_els; j++ ) {
      fprintf( outf, format, a[j]);
   }
   fprintf( outf, "\n");
}

void dWriteScalar( FILE *outf, char *text, double a, char *format)
{
   int j, k;

   fprintf( outf, "* %s\n", text);
   fprintf( outf, "S ");
   fprintf( outf, format, a);
   fprintf( outf, "\n");
}

void dvadd( double *a, int a_els, double *b, double *y)
{
   int j;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvadd");
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to 2nd vector: dvadd");
   if ( !valid_dvector_b( y ) ) 
      nrerror("Null pointer to result vector: dvadd");
#endif
   for ( j=1; j<=a_els; j++ ) {
      y[j] = a[j] + b[j];
   }
}   

void dvsub( double *a, int a_els, double *b, double *y)
{
   int j;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvsub");
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to 2nd vector: dvsub");
   if ( !valid_dvector_b( y ) ) 
      nrerror("Null pointer to result vector: dvsub");
#endif

   for ( j=1; j<=a_els; j++ ) {
      y[j] = a[j] - b[j];
   }
}   

double dvdot( double *a, int a_els, double *b)
{
   int j;
   double y;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvdot");
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to 2nd vector: dvdot");
#endif

   y = 0.0;
   for ( j=1; j<=a_els; j++ ) {
      y += a[j] * b[j];
   }
   return(y);
}   

#define TINY 1.0E-20
double dvmag( double *a, int a_els)
{
   int j;
   double y;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to vector: dvmag");
#endif

   y = 0.0;
   for ( j=1; j<=a_els; j++ ) {
      y += a[j] * a[j];
   }
   if ( y > TINY ) y = sqrt( y );
   return(y);
}   

double dvmnorm( double *a, int a_els)
{
   int j;
   double y;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to vector: dvmnorm");
#endif

   y = 0.0;
   for ( j=1; j<=a_els; j++ ) {
      y += a[j] * a[j];
   }
   if ( y > TINY ) {
      y = sqrt( y );
      for ( j=1; j<=a_els; j++ ) {
         a[j] = a[j]/y;
      }   
   }
   return(y);
}   
#undef TINY

void dvsmy( double *a, int a_els, double r, double *y)
{
   int j;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvsmy");
   if ( !valid_dvector_b( y ) ) 
      nrerror("Null pointer to result vector: dvsmy");
#endif

   for ( j=1; j<=a_els; j++ ) {
      y[j] = a[j] * r;
   }
}   

void dvpiv( double *a, int a_els, double r, double *b, double *y)
{
   int j;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvpiv");
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to 2nd vector: dvpiv");
   if ( !valid_dvector_b( y ) ) 
      nrerror("Null pointer to result vector: dvpiv");
#endif

   for ( j=1; j<=a_els; j++ ) {
      y[j] = a[j] * r + b[j];
   }
}   

/*--------------------------------------------------------------------
  Function for reading in a matrix or vector from a data file
  Format:
  Line 1:      M   nrow  {ncol}
  Line 2:      ncol real values
  ..
  ..
  Line 1+nrow: ncol real values

*/
#define BUF 1000


void dReadMatrix(FILE *datafile, double **y, int *rows, int *cols, int *error)
/* read a matrix (NR type) */
{
   char    line[BUF+2], *token;
   int     i, j, row, col, nrow, ncol, found, completed;
   double  r;
   char *fgr;

   *error = 0;
   nrow = 0;
   ncol = 0;
   completed = 0;
   found = 0;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( y ) ) 
      nrerror("Invalid matrix: dReadMatrix\n");   
#endif

   if (datafile == NULL) {
      fprintf(stderr,"ReadMatrix: data file not open\n");
      *error = 1;
      return;
   }
   row = 1;
   while ((fgr = fgets(line, BUF, datafile)) != NULL && !found) {
      if ( strncmp(line, "M", 1) == 0) {
         found = 1;                 
      }
      if ( found ) {
         /*          fprintf(stderr,"ReadMatrix: %s",line); */
			token = strtok(line, " ,\t\n");
			token = strtok(NULL, " ,\t\n"); /* skip the 'M' */
			if (token != NULL) {
				nrow = atoi(token);
				token = strtok(NULL, " ,\t\n");
				if (token != NULL) {
					ncol = atoi(token);
				}
				else {
					ncol = 1;                 /* presume column matrix */
				}
			}
			else {
				fprintf(stderr,"ReadMatrix: Could not read *rows* in data");
				*error = 2;
			}
		}

	}
	if ( !found ) {
		fprintf(stderr,"ReadMatrix: M for matrix not found!\n");
		*error = 3;
		exit(1);
	}

	while (fgr != NULL && !completed) {
		if (strlen(line) == 0 ||
			strncmp(line, "*", 1) == 0 ||
			strncmp(line, "\n", 1) == 0) {
			;    /* Blank line or comment */
		}
		else { /* read in the data */

			if (row <= nrow && *error == 0) {
				token = strtok(line, " ,\t\n");
				i = 1;
				while (token != NULL && i <= ncol) {
					sscanf(token, "%lf", &r);
					y[row][i] = r;
					/* printf("%s=y[%d][%d]: %lf\n",token,row,i,r);  */
					token = strtok(NULL, " ,\t\n");
					i++;
				}
				/* getchar(); */
				if ( row == nrow )  {
					*rows = nrow;
					*cols = ncol;
					return;
				}
				row++;
			}
		}     /* not a comment */
		fgr = fgets(line, BUF, datafile);
	}       /* while not EOF */
	if ( !completed ) {
		fprintf(stderr,"ReadMatrix: more data expected!\n");
		*error = 4;
		return;
	}
}

/*--------------------------------------------------------------------
  Function for reading in a vector from a data file
  Format:
  Line 1:      V   nels
  Line 2:      nels real values
*/
void dReadVector(FILE *datafile, double *y, int *els, int *error)
/* read a vector (NR type) */
{
	char    line[BUF+2], *token;
	int     i, j, nels, completed, found;
	double  r;
	char *fgr;

	*error = 0;
	nels = 0;
	completed = 0;
	found = 0;

	if (datafile == NULL) {
		fprintf(stderr,"ReadVector: data file not open\n");
		*error = 1;
		return;
	}
	while ((fgr = fgets(line, BUF, datafile)) != NULL && !found) {
		if ( strncmp(line, "V", 1) == 0) {
			found = 1;
		}
		if ( found ) {
			/*          fprintf(stderr,"ReadVector:%s",line); */
			token = strtok(line, " ,\t\n");
			token = strtok(NULL, " ,\t\n"); /* skip the 'V' */
			if (token != NULL) {
				nels = atoi(token);
			}
			else {
				fprintf(stderr,"ReadVector: Could not read *nels* in data");
				*error = 2;
				exit(1);
			}
		}
	}
	if ( !found ) {
		fprintf(stderr,"ReadVector: V for Vector not found!\n");
		*error = 3;
		exit(1);
	}

	while (fgr != NULL && !completed) {

		if (strlen(line) == 0 ||
			strncmp(line, "*", 1) == 0 ||
			strncmp(line, "\n", 1) == 0) {
			;    /* Blank line or comment */
		}
		else { /* read in the data */
			if (*error == 0) {
				token = strtok(line, " ,\t\n");
				i = 1;
				while (token != NULL && i <= nels) {
					sscanf(token, "%lf", &r);
					y[i] = r;
					/* printf("%s=y[%d]: %lf\n",token,i,r); */
					token = strtok(NULL, " ,\t\n");
					i++;
				}
				if ( i > nels ) {
					*els = nels;
					return;
				}
				/* getchar(); */
			}
		}      /* not a comment */
		fgr = fgets(line, BUF, datafile);
	}        /* while not EOF */
	if ( !completed ) {
		fprintf(stderr,"ReadVector: more data expected!\n");
		*error = 5;
		return;
	}
}

/*--------------------------------------------------------------------
  Function for reading in a scalar from a data file
  Format:
  Line 1:      S   value
*/
double dReadScalar(FILE *datafile, int *error)
{
	char    line[BUF+2], *token;
	int     i, found;
	double  r;

	*error = 0;
	found = 0;

	if (datafile == NULL) {
		fprintf(stderr,"ReadScalar: data file not open\n");
		*error = -1;
		return(0.0);
	}
	while (fgets(line, BUF, datafile) != NULL && !found) {
		if ( strncmp(line, "S", 1) == 0) {
			found = 1;
		}
		if ( found ) {
			/*          fprintf(stderr,"ReadScalar: %s",line); */
			token = strtok(line, " ,\t\n");
			token = strtok(NULL, " ,\t\n"); /* skip the 'S' */
			if (token != NULL) {
				r = atof(token);
			}
			else {
				fprintf(stderr,"ReadScalar: Could not read scalar in data");
				*error = -2;
				exit(1);
			}
		}
	}
	if ( !found ) {
		fprintf(stderr,"ReadVector: S for Scalar not found!\n");
		*error = 10;
		exit(1);
	}
	return(r);
}


void dmcopy( double **a, int a_rows, int a_cols, double **b) /* copy a matrix */
{
	int i, j;

#ifdef V_CHECK
	if ( !valid_dmatrix_b( a ) )
		nrerror("Invalid matrix: dmcopy\n");
	if ( !valid_dmatrix_b( b ) )
      nrerror("Invalid destination matrix: dmcopy\n");    
#endif

   for ( i=1; i<=a_rows; i++ )
      for ( j=1; j<=a_cols; j++ ) b[i][j] = a[i][j];
}

void dvcopy( double *a, int a_els, double *y) /* copy a vector */
{
   int i;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvcopy");
   if ( !valid_dvector_b( y ) ) 
      nrerror("Null pointer to 2nd vector: dvcopy");
#endif  

   for ( i=1; i<=a_els; i++ ) y[i] = a[i];
}

double dvcomp( double *a, int a_els, double *b) /* compare two vectors */
{
   int i;
   double ma, mb, da, sa, ra;

#ifdef V_CHECK
   if ( !valid_dvector_b( a ) ) 
      nrerror("Null pointer to 1st vector: dvcomp");
   if ( !valid_dvector_b( b ) ) 
      nrerror("Null pointer to 2nd vector: dvcomp");
#endif  

   ma = dvmag( a, a_els );
   mb = dvmag( b, a_els );
   da = 0.0;
   for (i=1; i<=a_els; i++) da += DSQR(a[i] - b[i]);   
   sa = sqrt( ma*mb );
   da = sqrt( da );
   if ( ma > mb ) ra = 1.0 - mb/ma;
   else ra = 1.0 - ma/mb;
   /*  printf("Mag rat %lf ", ra); */
   if ((da = da/sa) > ra ) ra = da;
   /*  printf("diff %lf ", da); */
   return(ra);
}

void dgetcolumn(double **G, int col, double *v, int nels) /* get a column from a matrix */
{
   int i;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( G ) ) 
      nrerror("Invalid matrix: dgetcolumn\n");    
   if ( !valid_dvector_b( v ) ) 
      nrerror("Null pointer to vector: dgetcolumn");
#endif

   for ( i=1; i<=nels; i++ ) v[i] = G[i][col];
}

void dputcolumn(double *v, int nels, double **G, int col) /* put a column into a matrix */
{
   int i;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( G ) ) 
      nrerror("Invalid matrix: dputcolumn\n");    
   if ( !valid_dvector_b( v ) ) 
      nrerror("Null pointer to vector: dputcolumn");
#endif

   for ( i=1; i<=nels; i++ ) G[i][col] = v[i];
}



#ifdef TEST
main()
{
   /* test program for above utility routines */

   double **a, **b, **c, **bT;
   double *x, *y, *z;
   FILE *infile, *outfile;
   int a_rows, a_cols, b_rows, b_cols, errors, xn, yn;

   infile = fopen("mat.in", "r");
   outfile = fopen("mat.dat", "w");

   a = dReadMatrix( infile, &a_rows, &a_cols, &errors);
   b = dReadMatrix( infile, &b_rows, &b_cols, &errors);
   x = dReadVector( infile, &xn, &errors);
   y = dReadVector( infile, &yn, &errors);
   getchar();

   dmdump( stdout, "Matrix A", a, a_rows, a_cols, "%8.2lf");
   dmdump( stdout, "Matrix B", b, b_rows, b_cols, "%8.2lf");
   dvdump( stdout, "Vector x", x, xn, "%8.2lf");
   dvdump( stdout, "Vector y", y, yn, "%8.2lf");
   z = dvector( 1, xn );
   dvadd( x, xn, y, z );
   dvdump( stdout, "x + y", z, xn, "%8.2lf");
   dvsub( x, xn, y, z );
   dvdump( stdout, "x - y", z, xn, "%8.2lf");
   dvsmy( x, xn, 2.0, z );
   dvdump( stdout, "2x", z, xn, "%8.2lf");
   printf("Magnitude of 2x: %7.2lf\n", dvmag( z, xn ));
   printf("dot product x.y: %7.2lf\n", dvdot( x, xn, y));

   dmvmult( a, a_rows, a_cols, x, xn, z );
   dvdump( stdout, "Ax", z, xn, "%8.2lf");

   c = dmatrix( 1, a_rows, 1, b_cols );    
   bT = dmatrix( 1, b_cols, 1, b_rows );
   dmtranspose( b, b_rows, b_cols, bT);
   dmdump( stdout, "Matrix B (transposed)", bT, b_cols, b_rows, "%8.2lf");
   dmmult( a, a_rows, a_cols, bT, b_cols, b_rows, c);
   dmdump( stdout, "Matrix AB", c, a_rows, b_rows, "%8.2lf");

   /*  dmfillUT( a, a_rows, a_cols );
       dmdump( stdout, "Symmetrified matrix A", a, a_rows, a_cols, "%8.2lf"); */

   free_dmatrix( a, 1, a_rows, 1, a_cols);
   free_dmatrix( b, 1, b_rows, 1, b_cols);
   free_dmatrix( c, 1, a_rows, 1, b_cols);
   free_dvector( x, 1, xn );
   free_dvector( y, 1, yn );
}
#endif
