/* nrlin.c ****************************************************************

  NR Library Package - Definitions and Linear Algebra Functions

  James Trevelyan,  University of Western Australia
  Revision 2   January 1996

**************************************************************************/


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#ifdef __TURBOC__
#include <conio.h>
#endif

#include "nrlin.h"

/* Validity checking - use V_CHECK */
#define V_CHECK

/* Keep track of memory used - only for dmatrix and dvector calls at the moment */

static long mem_t = 0;
static int mats = 0;
static int dmats = 0;
static int vecs = 0;
static int dvecs = 0;

void using_mem( long space )
{
	mem_t += space;
}

long mem_used( void )
{
	return (mem_t);
}

void reportmemory( FILE *outfile )
{
	 fprintf( outfile, "Memory: %ld bytes\n", mem_t);
	 fprintf( outfile, "%d dmatrices, %d dvectors, %d matrices, %d vectors\n",
				 dmats, dvecs, mats, vecs );
}
void choldc(float **a, int n, float p[])
{
	nrerror("Obtain \"choldc\" source from NR diskette");
}


void cholsl(float **a, int n, float p[], float b[], float x[])
{
	nrerror("Obtain \"cholsl\" source from NR diskette");
}


void lubksb(float **a, int n, int *indx, float b[])
{
	nrerror("Obtain \"lubksb\" source from NR diskette");
}


void ludcmp(float **a, int n, int *indx, float *d)
{
	nrerror("Obtain \"ludcmp\" source from NR diskette");
}

void dcholdc(double **a, int n, double p[])
{
	nrerror("Obtain \"choldc\" source from NR diskette");
}


void dcholsl(double **a, int n, double p[], double b[], double x[])
{
	nrerror("Obtain \"cholsl\" source from NR diskette");
}


void dlubksb(double **a, int n, int *indx, double b[])
{
	nrerror("Obtain \"lubksb\" source from NR diskette");
}


void dludcmp(double **a, int n, int *indx, double *d)
{
	nrerror("Obtain \"ludcmp\" source from NR diskette");
}




#define NR_END 1
#define NR_TEST 1
#define FREE_ARG char*
static int alt_handler_defined = 0;
static void (*alt_error_handler)(char error_text[]);

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	if (alt_handler_defined) {
      (*alt_error_handler)(error_text);
   } 
   else {
      fprintf(stderr,"Numerical Recipes run-time error...\n");
      fprintf(stderr,"%s\n",error_text);
      fprintf(stderr,"Press <ENTER> to return to system...");
      getchar();
      exit(1);
   }
}

void nrerror_handler( void(*handler)(char error_text[]) )
{
   alt_handler_defined = 1;
   alt_error_handler = handler;
}

float *vector(I_ARG_T nl, I_ARG_T nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END+NR_TEST)*sizeof(float)));
	using_mem( (nh-nl+1+NR_END+NR_TEST)*sizeof(float) );
	if (!v) nrerror("allocation failure in dvector()");
	v -= nl;
	v += NR_END;
	v[nl-1] = -322.0;
	v[nh+1] = -722.0;
	vecs++;
	return (v);
}

int *ivector(I_ARG_T nl, I_ARG_T nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
   int *v;

   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(I_ARG_T nl, I_ARG_T nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

long *lvector(I_ARG_T nl, I_ARG_T nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	long *v;

	v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

void flag_dvector( double *v, I_ARG_T nh)
{
	v[0] = -322.0;
	v[nh+1] = -722.0;
}

double *dvector(I_ARG_T nl, I_ARG_T nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END+NR_TEST)*sizeof(double)));
	using_mem( (nh-nl+1+NR_END+NR_TEST)*sizeof(double) );
	if (!v) nrerror("allocation failure in dvector()");
	v -= nl;
	v += NR_END;
	v[nl-1] = -322.0;
	v[nh+1] = -722.0;
	dvecs++;
	return (v);
}

float **matrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	I_ARG_T i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END+1)*sizeof(float*)));
	using_mem( (nrow+NR_END+1)*sizeof(float*) );
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END + 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END+NR_TEST)*sizeof(float)));
	using_mem( (nrow*ncol+NR_END+NR_TEST)*sizeof(float) );
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	m[nrl-2] = (float *)0x555;
	m[nrl-1] = m[nrl];  /* covers m[0][0] mistake */
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	m[nrl][0]=-422.0;
	m[nrh][nch+1]=-822.0;

	/* return pointer to array of pointers to rows */
	mats++;
	return m;
}

double **dmatrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	I_ARG_T i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END+1)*sizeof(double*)));
	using_mem( (nrow+NR_END+1)*sizeof(double*) );
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END + 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END+NR_TEST)*sizeof(double)));
	using_mem( (nrow*ncol+NR_END+NR_TEST)*sizeof(double) );
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	m[nrl-2] = (double *)0x555;
	m[nrl-1] = m[nrl];   /* covers m[0][0] mistake */
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	m[nrl][0]=-422.0;
	m[nrh][nch+1]=-822.0;
	/* return pointer to array of pointers to rows */
	dmats++;
	return m;
}

int **imatrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	I_ARG_T i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END+NR_TEST)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
   return m;
}

float **submatrix(float **a, I_ARG_T oldrl, I_ARG_T oldrh, I_ARG_T oldcl, I_ARG_T oldch,
    I_ARG_T newrl, I_ARG_T newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
   I_ARG_T i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
   float **m;

   /* allocate array of pointers to rows */
   m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
   if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
   m -= newrl;

   /* set pointers to rows */
   for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}

float **convert_matrix(float *a, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	I_ARG_T i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
   float **m;

   /* allocate pointers to rows */
   m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
   if (!m) nrerror("allocation failure in convert_matrix()");
   m += NR_END;
   m -= nrl;

   /* set pointers to rows */
   m[nrl]=a-ncl;
   for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
   /* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch, I_ARG_T ndl, I_ARG_T ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
   I_ARG_T i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   float ***t;

   /* allocate pointers to pointers to rows */
   t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
   t += NR_END;
   t -= nrl;

   /* allocate pointers to rows and set pointers to them */
   t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
   if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
   t[nrl] += NR_END;
	t[nrl] -= ncl;

   /* allocate rows and set pointers to them */
   t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
   if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++) {
      t[i]=t[i-1]+ncol;
      t[i][ncl]=t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, I_ARG_T nl, I_ARG_T nh)
/* free a float vector allocated with vector() */
{
	if ( valid_vector( v, nl, nh ) ) {
		free((FREE_ARG) (v+nl-NR_END));
		using_mem( -(nh-nl+1+NR_END+NR_TEST)*sizeof(float) );
		vecs--;
	}
	else
		nrerror("Invalid vector pointer: free_vector");
}

void free_ivector(int *v, I_ARG_T nl, I_ARG_T nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, I_ARG_T nl, I_ARG_T nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(long *v, I_ARG_T nl, I_ARG_T nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, I_ARG_T nl, I_ARG_T nh)
/* free a double vector allocated with dvector() */
{
	if ( valid_dvector( v, nl, nh ) ) {
		free((FREE_ARG) (v+nl-NR_END));
		using_mem( -(long)(nh-nl+1+NR_END+NR_TEST)*sizeof(double) );
		dvecs--;
	}
	else
		nrerror("Invalid vector pointer: free_vector");
}

void free_matrix(float **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* free a float matrix allocated by matrix() */
{
	I_ARG_T i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	if ( valid_matrix( m, nrl, nrh, ncl, nch ) ) {
		m[nrl][ncl-1] = 0.0;
		m[nrh][nch+1] = 0.0;
		*(m-1) = (float *)NULL;
		free((FREE_ARG) (m[nrl]+ncl-NR_END));
		free((FREE_ARG) (m+nrl-NR_END-1));
		using_mem(-(long)(nrow+NR_END+1)*sizeof(float*) );
		using_mem(-(long)(nrow*ncol+NR_END+NR_TEST)*sizeof(float) );
		mats--;
	}
	else {
		nrerror("Invalid pointer to matrix: free_matrix");
	}
}

void free_dmatrix(double **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* free a double matrix allocated by dmatrix() */
{
	I_ARG_T i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	if ( valid_dmatrix( m, nrl, nrh, ncl, nch ) ) {
		m[nrl][ncl-1] = 0.0;
		m[nrh][nch+1] = 0.0;
		*(m-1) = (double *)NULL;
		free((FREE_ARG) (m[nrl]+ncl-NR_END));
		free((FREE_ARG) (m+nrl-NR_END-1));
		using_mem(-(long)(nrow+NR_END+1)*sizeof(double*) );
		using_mem(-(long)(nrow*ncol+NR_END+NR_TEST)*sizeof(double) );
		dmats--;
	}
	else {
		nrerror("Invalid pointer to dmatrix: free_dmatrix");
	}
}

void free_imatrix(int **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch,
	 I_ARG_T ndl, I_ARG_T sndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

/* matrix inversion ref - Page 48 */


void dinverse( double **a, int n, double **y )
/*  Find inverse of 'a' (decomposed in process!) and return as 'y' */
{
	double d, *col;
	int i, j, *indx /* integer vector */;

#ifdef V_CHECK
	if ( !valid_dmatrix_b( a ) )
		nrerror("Invalid input matrix: dinverse");
	if ( !valid_dmatrix_b( y ) )
		nrerror("Invalid output matrix: dinverse");
#endif
	indx = ivector( 1, n );
	col = dvector( 1, n );
	dludcmp( a, n, indx, &d);
	for ( j=1; j<=n; j++ ) {
		for ( i=1; i<=n; i++ ) col[i] = 0.0;
		col[j] = 1.0;
		dlubksb( a, n, indx, col );
		for ( i=1; i<=n; i++ ) y[i][j] = col[i];
	}
	free_ivector( indx, 1, n );
   free_dvector( col, 1, n );
}

/* if intending to compute inverse(A) * B, use columns of B in 'col' above instead
   of unit vectors as shown - more accurate and faster */

void dinverse_mult( double **a, int a_rows, double **b, int b_cols, double **y )
/*  Find inverse of 'a' (decomposed in process!) times 'b' and return result as 'y' */
{
   double d, *col;
   int i, j, *indx /* integer vector */;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) )
      nrerror("Invalid input matrix: dinverse_mult");
   if ( !valid_dmatrix_b( b ) )
      nrerror("Invalid output matrix: dinverse_mult");
   if ( !valid_dmatrix_b( y ) )
      nrerror("Invalid output matrix: dinverse_mult");
#endif
   indx = ivector( 1, a_rows );
   col = dvector( 1, a_rows );
   dludcmp( a, a_rows, indx, &d);  
   for ( j=1; j<=b_cols; j++ ) {
      for ( i=1; i<=a_rows; i++ ) col[i] = b[i][j];
      dlubksb( a, a_rows, indx, col );
      for ( i=1; i<=a_rows; i++ ) y[i][j] = col[i];
   }
   free_ivector( indx, 1, a_rows );
   free_dvector( col, 1, a_rows );
}




/* positive definite symmetric matrix inversion ref - Page 97, 98 */


void dPDSinverse( double **a, int n, double **y )
/*  Find inverse of 'a' (decomposed in process!) and return as 'y' */
{
   double sum, *p, *col, *yr;
   int i, j;

#ifdef V_CHECK
   if ( !valid_dmatrix_b( a ) )
      nrerror("Invalid input matrix: dPDSinverse");
   if ( !valid_dmatrix_b( y ) )
      nrerror("Invalid output matrix: dPDSinverse");
#endif
   p = dvector( 1, n );
   col = dvector( 1, n );
   yr = dvector( 1, n );
   dcholdc( a, n, p);  
   for ( j=1; j<=n; j++ ) {
      for ( i=1; i<=n; i++ ) col[i] = 0.0;
      col[j] = 1.0;
      dcholsl( a, n, p, col, yr );
      for ( i=1; i<=n; i++ ) y[i][j] = yr[i];
   }
   free_dvector( yr, 1, n );
   free_dvector( col, 1, n );
   free_dvector( p, 1, n );
}


void dPDS_L_inverse( double **a, int n, double **y )
/*  Find inverse of L (decomposed a) and return as 'y' */
{
   double sum, *p;
   int i, j, k;
#ifdef V_CHECK
   if ( !valid_dmatrix_b( a) )
      nrerror("Invalid input matrix: dPDS_L_inverse");
   if ( !valid_dmatrix_b( y ) )
      nrerror("Invalid output matrix: dPDS_L_inverse");
#endif
   p = dvector( 1, n );
   dcholdc( a, n, p);  
   for ( i=1; i<=n; i++ ) {
      a[i][i] = 1.0/p[i];
      for ( j=i+1; j<=n; j++ ) {
         sum = 0.0;
         for (k=i; k<j; k++) sum -= a[j][k]*a[k][i];
         a[j][i] = sum/p[j];
      }
   }
   free_dvector( p, 1, n );
}
