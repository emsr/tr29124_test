/* NRLIN.H - include file  */

#ifndef _NRLIN_H_
#define _NRLIN_H_

#define I_ARG_T int

void choldc(float **a, int n, float p[]);
void cholsl(float **a, int n, float p[], float b[], float x[]);
void lubksb(float **a, int n, int *indx, float b[]);
void ludcmp(float **a, int n, int *indx, float *d);

void dcholdc(double **a, int n, double p[]);
void dcholsl(double **a, int n, double p[], double b[], double x[]);
void dlubksb(double **a, int n, int *indx, double b[]);
void dludcmp(double **a, int n, int *indx, double *d);


#define validate_dvector(a,b,c) a[(b)-1]=-322.0;a[(c)+1]=-722.0
#define validate_vector(a,b,c) a[(b)-1]=-322.0;a[(c)+1]=-722.0
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


static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : (sqrarg)*(sqrarg))

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : (dsqrarg)*(dsqrarg))

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
		  (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
		  (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

long mem_used(void);
void reportmemory( FILE *outfile );
void using_mem( long space );
void nrerror(char error_text[]);
float *vector(I_ARG_T nl, I_ARG_T nh);
int *ivector(I_ARG_T nl, I_ARG_T nh);
unsigned char *cvector(I_ARG_T nl, I_ARG_T nh);
long *lvector(I_ARG_T nl, I_ARG_T nh);
double *dvector(I_ARG_T nl, I_ARG_T nh);
void flag_dvector( double *v, I_ARG_T nh);
float **matrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
double **dmatrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
int **imatrix(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
float **submatrix(float **a, I_ARG_T oldrl, I_ARG_T oldrh, I_ARG_T oldcl, I_ARG_T oldch,
    I_ARG_T newrl, I_ARG_T newcl);
float **convert_matrix(float *a, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
float ***f3tensor(I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch, I_ARG_T ndl, I_ARG_T ndh);
void free_vector(float *v, I_ARG_T nl, I_ARG_T nh);
void free_ivector(int *v, I_ARG_T nl, I_ARG_T nh);
void free_cvector(unsigned char *v, I_ARG_T nl, I_ARG_T nh);
void free_lvector(long *v, I_ARG_T nl, I_ARG_T nh);
void free_dvector(double *v, I_ARG_T nl, I_ARG_T nh);
void free_matrix(float **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
void free_dmatrix(double **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
void free_imatrix(int **m, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
void free_submatrix(float **b, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
void free_convert_matrix(float **b, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch);
void free_f3tensor(float ***t, I_ARG_T nrl, I_ARG_T nrh, I_ARG_T ncl, I_ARG_T nch,
    I_ARG_T ndl, I_ARG_T ndh);

void dinverse( double **a, int n, double **y );
void dinverse_mult( double **a, int a_rows, double **b, int b_cols, double **y );
void dPDSinverse( double **a, int n, double **y );
void dPDS_L_inverse( double **a, int n, double **y );


/* sources in UTIL.C */

void dmmult( double **a, int a_rows, int a_cols, 
             double **b, int b_rows, int b_cols, double **y);
void dmvmult( double **a, int a_rows, int a_cols, double *b, int b_els, double *y);
void dmadd( double **a, int a_rows, int a_cols, double **b, double **y);
void dmsub( double **a, int a_rows, int a_cols, double **b, double **y);
void dmsmy( double **a, int a_rows, int a_cols, double r, double **y);
void dmtranspose( double **a, int a_rows, int a_cols, double **y);
void dmfillUT( double **a, int a_rows, int a_cols);

void dmdump( FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format);
void dvdump( FILE *outf, char *text, double *a, int a_els, char *format);
void dvadd( double *a, int a_els, double *b, double *y);
void dvsub( double *a, int a_els, double *b, double *y);
double dvdot( double *a, int a_els, double *b);
double dvmag( double *a, int a_els);
double dvmnorm( double *a, int a_els);
void dvsmy( double *a, int a_els, double r, double *y);
void dvpiv( double *a, int a_els, double r, double *b, double *y);
void dReadMatrix(FILE *datafile, double **y, int *rows, int *cols, int *error);
void dReadVector(FILE *datafile, double *y, int *els, int *error);
double dReadScalar(FILE *datafile, int *error);
void dWriteMatrix( FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format);
void dWriteVector( FILE *outf, char *text, double *a, int a_els, char *format);
void dWriteScalar( FILE *outf, char *text, double a, char *format);

void dmcopy( double **a, int a_rows, int a_cols, double **b); /* copy a matrix */
void dvcopy( double *a, int a_els, double *b);                /* copy a vector */
double dvcomp( double *a, int a_els, double *b);              /* compare two vectors */
void dgetcolumn(double **G, int col, double *v, int nels);    /* get a column from a matrix */
void dputcolumn(double *v, int nels, double **G, int col);    /* put a column into a matrix */

/* frame operations  - FRAME.C */
double **dframe( void );
double *d4vector( void );
#define free_d4vector( y )  free_dvector( y, 1, 4 )
#define free_dframe( y )    free_dmatrix( y, 1, 4, 1, 4 )
#define valid_d4vector( y ) valid_dvector( y, 1, 4 )
#define valid_d4vector_b( y ) valid_dvector_b( y )
#define valid_dframe( y )   valid_dmatrix( y, 1, 4, 1, 4 )

void chainmult( double **a, double **b );
void rotframe( double **a, int axis, double theta); /* make a rotation transform */
void transframe( double **a, double dx, double dy, double dz); 
void orthog( double **a );                          /* orthogonalize a frame */
void orthogKI( double **a );                          /* orthogonalize a frame */
void dframeinverse( double **E, double **Ein ); /* form inverse frame */
void dcross4 ( double *a, double *b, double *c );   /* cross two vectors */
double column_dot( double **G1, int col1, double **G2, int col2 );
void make_d_omega( double **G1, double **G2, double *w);
void compare_link_frames( double **G1, double **G2, double *wp, double *wr, double *ww );
void compare_frames( double **G1, double **G2, double *wp, double *wr, double *ww );
void FindAxis(double **T1, double **T2, double *axis, 
                double *sphi, double *cphi,int *twitching);
void RotateVector(double *v, double *axis, double sphi, double cphi, double *y);
void RotateFrame(double **T1, double *axis, double sphi, double cphi, double **T2);

#ifdef NR_CHECK
/* include definitions for wrapper functions which check parameters
	and locate error calls */
#define free_dmatrix(a,b,c,d,e) \
		free_dmatrix_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define free_dvector(a,b,c) \
		free_dvector_DNR(__LINE__,__FILE__,a,b,c)

#define dmatrix(a,b,c,d) \
		dmatrix_DNR(__LINE__,__FILE__,a,b,c,d)
#define dvector(a,b) \
		dvector_DNR(__LINE__,__FILE__,a,b)

#define dmdump(a,b,c,d,e,f) \
		dmdump_DNR(__LINE__, __FILE__,a,b,c,d,e,f)
#define dvdump(a,b,c,d,e) \
		dvdump_DNR(__LINE__, __FILE__,a,b,c,d,e)
#define dWriteMatrix(a,b,c,d,e,f) \
		dWriteMatrix_DNR(__LINE__, __FILE__,a,b,c,d,e,f)
#define dWriteVector(a,b,c,d,e) \
		dWriteVector_DNR(__LINE__, __FILE__,a,b,c,d,e)
#define dWriteScalar(a,b,c,d) \
		dWriteScalar_DNR(__LINE__, __FILE__,a,b,c,d)
#define dReadMatrix(a,b,c,d,e) \
		dReadMatrix_DNR(__LINE__, __FILE__,a,b,c,d,e)
#define dReadVector(a,b,c,d) \
		dReadVector_DNR(__LINE__, __FILE__,a,b,c,d)
#define dReadScalar(a,b) \
		dReadScalar_DNR(__LINE__, __FILE__,a,b)

#define dmmult(a,b,c,d,e,f,g)\
		dmmult_DNR(__LINE__,__FILE__,a,b,c,d,e,f,g)
#define dmvmult(a,b,c,d,e,f) \
		dmvmult_DNR(__LINE__,__FILE__,a,b,c,d,e,f)
#define dmadd(a,b,c,d,e) \
		dmadd_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dmsub(a,b,c,d,e) \
		dmsub_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dmsmy(a,b,c,d,e) \
		 dmsmy_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dmtranspose(a,b,c,d) \
		dmtranspose_DNR(__LINE__,__FILE__,a,b,c,d)
#define dmfullUT(a,b,c) \
		dmfillUT_DNR(__LINE__,__FILE__,a,b,c)

#define dvadd(a,b,c,d) \
		 dvadd_DNR(__LINE__,__FILE__,a,b,c,d)
#define dvsub(a,b,c,d) \
		 dvsub_DNR(__LINE__,__FILE__,a,b,c,d)
#define dvsmy(a,b,c,d) \
		 dvsmy_DNR(__LINE__,__FILE__,a,b,c,d)
#define dvpiv(a,b,c,d,e) \
		 dvpiv_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dvdot(a,b,c) \
		 dvdot_DNR(__LINE__,__FILE__,a,b,c)
#define dvmag(a,b) \
		 dvmag_DNR(__LINE__,__FILE__,a,b)
#define dvmnorm(a,b) \
		 dvmnorm_DNR(__LINE__,__FILE__,a,b)

#define dmcopy(a,b,c,d) \
		 dmcopy_DNR(__LINE__,__FILE__,a,b,c,d)
#define dvcopy(a,b,c) \
		 dvcopy_DNR(__LINE__,__FILE__,a,b,c)
#define dvcomp(a,b,c) \
		 dvcomp_DNR(__LINE__,__FILE__,a,b,c)
#define dgetcolumn(a,b,c,d) \
		 dgetcolumn_DNR(__LINE__,__FILE__,a,b,c,d)
#define dputcolumn(a,b,c,d) \
		 dputcolumn_DNR(__LINE__,__FILE__,a,b,c,d)

void free_dmatrix_DNR( int linref, char *fileref,
					 double **a, int nrl, int nrh, int ncl, int nch);
void free_dvector_DNR( int linref, char *fileref,
					 double *a, int nvl, int nvh);

double **dmatrix_DNR( int linref, char *fileref,
					 int nrl, int nrh, int ncl, int nch );
double *dvector_DNR( int linref, char *fileref,
					 int ncl, int nch );

void dmdump_DNR( int linref, char *fileref,
					FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format);
void dvdump_DNR( int linref, char *fileref,
					FILE *outf, char *text, double *a, int a_els, char *format);
void dWriteMatrix_DNR( int linref, char *fileref,
					FILE *outf, char *text, double **a, int a_rows, int a_cols, char *format);
void dWriteVector_DNR( int linref, char *fileref,
					FILE *outf, char *text, double *a, int a_els, char *format);
void dWriteScalar_DNR( int linref, char *fileref,
					FILE *outf, char *text, double a, char *format);
void dReadMatrix_DNR( int linref, char *fileref,
				  FILE *datafile, double **y, int *rows, int *cols, int *error);
void dReadVector_DNR( int linref, char *fileref,
				  FILE *datafile, double *y, int *els, int *error);
double dReadScalar_DNR( int linref, char *fileref,
				  FILE *datafile, int *error);

void dmmult_DNR( int linref, char *fileref,
				 double **a, int a_rows, int a_cols,
				 double **b, int b_rows, int b_cols, double **y);
void dmvmult_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double *b, int b_els, double *y);
void dmadd_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double **b, double **y);
void dmsub_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double **b, double **y);
void dmsmy_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double r, double **y);
void dmtranspose_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double **y);
void dmfillUT_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols);

void dvadd_DNR( int linref, char *fileref,
					double *a, int a_els, double *b, double *y);
void dvsub_DNR( int linref, char *fileref,
					double *a, int a_els, double *b, double *y);
double dvdot_DNR( int linref, char *fileref,
					double *a, int a_els, double *b);
double dvmag_DNR( int linref, char *fileref,
					double *a, int a_els);
double dvmnorm_DNR( int linref, char *fileref,
					double *a, int a_els);
void dvsmy_DNR( int linref, char *fileref,
					double *a, int a_els, double r, double *y);
void dvpiv_DNR( int linref, char *fileref,
					double *a, int a_els, double r, double *b, double *y);

void dmcopy_DNR( int linref, char *fileref,
					double **a, int a_rows, int a_cols, double **b); /* copy a matrix */
void dvcopy_DNR( int linref, char *fileref,
					double *a, int a_els, double *b);                /* copy a vector */
double dvcomp_DNR( int linref, char *fileref,
					double *a, int a_els, double *b);              /* compare two vectors */
void dgetcolumn_DNR( int linref, char *fileref,
				  double **G, int col, double *v, int nels);    /* get a column from a matrix */
void dputcolumn_DNR( int linref, char *fileref,
				  double *v, int nels, double **G, int col);    /* put a column into a matrix */

void dinverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y );
void dinverse_mult_DNR( int lineref, char *fileref,
					double **a, int a_rows, double **b, int b_cols, double **y );
void dPDSinverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y );
void dPDS_L_inverse_DNR( int lineref, char *fileref,
					double **a, int n, double **y );
void choldc_DNR( int lineref, char *fileref,
				  float **a, int n, float p[]);
void cholsl_DNR( int lineref, char *fileref,
				  float **a, int n, float p[], float b[], float x[]);
void lubksb_DNR( int lineref, char *fileref,
				  float **a, int n, int *indx, float b[]);
void ludcmp_DNR( int lineref, char *fileref,
				  float **a, int n, int *indx, float *d);

void dcholdc_DNR( int lineref, char *fileref,
				  double **a, int n, double p[]);
void dcholsl_DNR( int lineref, char *fileref,
				  double **a, int n, double p[], double b[], double x[]);
void dlubksb_DNR( int lineref, char *fileref,
				  double **a, int n, int *indx, double b[]);
void dludcmp_DNR( int lineref, char *fileref,
				  double **a, int n, int *indx, double *d);
#define dinverse(a,b,c) \
		  dinverse_DNR(__LINE__,__FILE__,a,b,c)
#define dinverse_mult(a,b,c,d,e) \
		  dinverse_mult_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dPDSinverse(a,b,c) \
		  dPDSinverse_DNR(__LINE__,__FILE__,a,b,c)
#define dPDS_L_inverse(a,b,c) \
		  dPDS_L_inverse_DNR(__LINE__,__FILE__,a,b,c)
#define choldc(a,b,c) \
		  choldc_DNR(__LINE__,__FILE__,a,b,c)
#define cholsl(a,b,c,d,e) \
		  cholsl_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define lubksb(a,b,c,d) \
		  lubksb_DNR(__LINE__,__FILE__,a,b,c,d)
#define ludcmp(a,b,c,d) \
		  ludcmp_DNR(__LINE__,__FILE__,a,b,c,d)
#define dcholdc(a,b,c) \
		  dcholdc_DNR(__LINE__,__FILE__,a,b,c)
#define dcholsl(a,b,c,d,e) \
		  dcholsl_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define dlubksb(a,b,c,d) \
		  dlubksb_DNR(__LINE__,__FILE__,a,b,c,d)
#define dludcmp(a,b,c,d) \
		  dludcmp_DNR(__LINE__,__FILE__,a,b,c,d)

double **dframe_DNR( int lineref, char *fileref );
double *d4vector_DNR( int lineref, char *fileref );

#define free_d4vector( y )  free_dvector( y, 1, 4 )
#define free_dframe( y )    free_dmatrix( y, 1, 4, 1, 4 )
#define valid_d4vector( y ) valid_dvector( y, 1, 4 )
#define valid_d4vector_b( y ) valid_dvector_b( y )
#define valid_dframe( y )   valid_dmatrix( y, 1, 4, 1, 4 )

void chainmult_DNR( int lineref, char *fileref,
					double **a, double **b );
void rotframe_DNR( int lineref, char *fileref,
					double **a, int axis, double theta); /* make a rotation transform */
void transframe_DNR( int lineref, char *fileref,
					double **a, double dx, double dy, double dz);
void orthog_DNR( int lineref, char *fileref,
					double **a );                          /* orthogonalize a frame */
void orthogKI_DNR( int lineref, char *fileref,
					double **a );                          /* orthogonalize a frame */
void dframeinverse_DNR( int lineref, char *fileref,
					double **E, double **Ein ); /* form inverse frame */
void dcross4_DNR( int lineref, char *fileref,
					double *a, double *b, double *c );   /* cross two vectors */
#define dframe() \
		  dframe_DNR(__LINE__,__FILE__)
#define d4vector() \
		  d4vector_DNR(__LINE__,__FILE__)
#define chainmult(a,b) \
		  chainmult_DNR(__LINE__,__FILE__,a,b)
#define rotframe(a,b,c) \
		  rotframe_DNR(__LINE__,__FILE__,a,b,c)
#define transframe(a,b,c,d) \
		  transframe_DNR(__LINE__,__FILE__,a,b,c,d)
#define othog(a) \
		  orthog_DNR(__LINE__,__FILE__,a)
#define orthogKI(a) \
		  orthogKI_DNR(__LINE__,__FILE__,a)
#define dframeinverse(a,b) \
		  dframeinverse_DNR(__LINE__,__FILE__,a,b)
#define dcross4(a,b,c) \
		  dcross4_DNR(__LINE__,__FILE__,a,b,c)
double column_dot_DNR( int lineref, char *fileref,
					double **G1, int col1, double **G2, int col2 );
void make_d_omega_DNR( int lineref, char *fileref,
					double **G1, double **G2, double *w);
void compare_link_frames_DNR( int lineref, char *fileref,
					double **G1, double **G2, double *wp, double *wr, double *ww );
void compare_frames_DNR( int lineref, char *fileref,
				  double **G1, double **G2, double *wp, double *wr, double *ww );
void FindAxis_DNR( int lineref, char *fileref,
				  double **T1, double **T2, double *axis,
					 double *sphi, double *cphi,int *twitching);
void RotateVector_DNR( int lineref, char *fileref,
				  double *v, double *axis, double sphi, double cphi, double *y);
void RotateFrame_DNR( int lineref, char *fileref,
				  double **T1, double *axis, double sphi, double cphi, double **T2);

#define dcolumn_dot(a,b,c,d) \
		  dcolumn_dot_DNR(__LINE__,__FILE__,a,b,c,d)
#define make_d_omega(a,b,c) \
		  make_d_omega_DNR(__LINE__,__FILE__,a,b,c)
#define compare_link_frames(a,b,c,d,e) \
		  compare_link_frames_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define compare_frames(a,b,c,d,e) \
		  compare_frames_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define FindAxis(a,b,c,d,e,f) \
		  FindAxis_DNR(__LINE__,__FILE__,a,b,c,d,e,f)
#define RotateVector(a,b,c,d,e) \
		  RotateVector_DNR(__LINE__,__FILE__,a,b,c,d,e)
#define RotateFrame(a,b,c,d,e) \
		  RotateFrame_DNR(__LINE__,__FILE__,a,b,c,d,e)
#endif /* NR_CHECK */
#endif /* _NR_H_ */

