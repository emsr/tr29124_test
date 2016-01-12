
#ifndef  __NRIC_H__
#define  __NRIC_H__

#ifdef __cplusplus
extern "C" {
#endif


#define PI2       6.2831853071795862
#define PI        3.1415926535897931
#define PIO2      1.5707963267948966
#define PIO3      1.0471975511965976
#define PIO4      0.7853981633974483

#define SQRT2     1.4142135623730950
#define SQRT2O2   0.7071067811865475
#define SQRT3     1.7320508075688773
#define SQRT3O2   0.8660254037844386

#define SQRTPI    1.7724538509055159
#define SQRTPIO2  1.2533141373155001
#define SQRTPI2   2.5066282746310005

#define LOG2      0.6931471805599453
#define LOGPI     1.1447298858494002

#define E         2.7182818284590452

#define GOLD      0.6180339887498949

#define GAMMA     0.5772156649015329



/*
 *    Complex in double precision.
 */
typedef struct FCOMPLEX {
    float r;
    float i;
} fcomplex;

/*
 *    Complex in single precision.
 */
typedef struct DCOMPLEX {
    double r;
    double i;
} dcomplex;



/*
 *    Three-dimensional vector in double precision.
 */
typedef struct DVECTOR3 {
    double x;
    double y;
    double z;
} dvector3;

/*
 *    Three-dimensional two-tensor in double precision.
 */
typedef struct DMATRIX3 {
    dvector3 x;
    dvector3 y;
    dvector3 z;
} dmatrix3;


/*
 *    Three-dimensional vector in single precision.
 */
typedef struct FVECTOR3 {
    float x;
    float y;
    float z;
} fvector3;

/*
 *    Three-dimensional two-tensor in single precision.
 */
typedef struct FMATRIX3 {
    fvector3 x;
    fvector3 y;
    fvector3 z;
} fmatrix3;

/*
 *    Two-dimensional vector in double precision.
 */
typedef struct DVECTOR2 {
    double x;
    double y;
} dvector2;

/*
 *    Three-dimensional two-tensor in double precision.
 */
typedef struct DMATRIX2 {
    dvector2 x;
    dvector2 y;
} dmatrix2;

/*
 *    Two dimensional vector in single precision.
 */
typedef struct FVECTOR2 {
    float x;
    float y;
} fvector2;

/*
 *    Two dimensional two-tensor in single precision.
 */
typedef struct FMATRIX2 {
    fvector2 x;
    fvector2 y;
} fmatrix2;


/*
 *    Three by two matrix in double precision.
 */
typedef struct DMATRIX32 {
    dvector2 x;
    dvector2 y;
    dvector2 z;
} dmatrix32;


/*
 *    Two by three matrix in double precision.
 */
typedef struct DMATRIX23 {
    dvector3 x;
    dvector3 y;
} dmatrix23;


/*
 *    Three by two matrix in single precision.
 */
typedef struct FMATRIX32 {
    fvector2 x;
    fvector2 y;
    fvector2 z;
} fmatrix32;


/*
 *    Two by three matrix in single precision.
 */
typedef struct FMATRIX23 {
    fvector3 x;
    fvector3 y;
} fmatrix23;


/*
const dvector3 DV3_ZERO = {0.0, 0.0, 0.0};
const fvector3 FV3_ZERO = {0.0, 0.0, 0.0};

const dvector2 DV2_ZERO = {0.0, 0.0};
const fvector2 FV2_ZERO = {0.0, 0.0};
*/


/*
 *    nric.c
 *
 *    Numerical Recipes error and warning routines.
 */
void    nrerror( const char *error_text);
void    nrwarning( const char *warning_text);


/*
 *    plot.c
 *
 *    Crude function plotting.
 */
void plot_func( double (*f)(double), double x1, double x2,
            const char *t1, const char *t2, const char *tx, const char *ty );


/*
 *    search_table.c
 *
 *    Searches an ordered table of numbers for the largest entry smaller than
 *    an input number.    Used in interpolation.
 */
void locate( double *xx, int n, double x, int *j );
void hunt( double *xx, int n, double x, int *j );



/*
 *    memory.c
 *
 *    Allocates memory for vectors and matrices and frees such vectors and matrices.
 */
unsigned char *ucvector( int nl, int nh);
void          free_ucvector( unsigned char *v, int nl, int nh);

char    *cvector( int nl, int nh);
void    free_cvector( char *v, int nl, int nh);
char    **cmatrix( int nrl, int nrh, int ncl, int nch);
void    free_cmatrix( char **m, int nrl, int nrh, int ncl, int nch);

unsigned int  *uivector( int nl, int nh);
void          free_uivector( unsigned int *v, int nl, int nh);

int     *ivector( int nl, int nh);
void    free_ivector( int *v, int nl, int nh);

unsigned long *ulvector( int nl, int nh);
void          free_ulvector( unsigned long *v, int nl, int nh);

long *lvector( int nl, int nh);
void free_lvector( long *v, int nl, int nh);

int     **imatrix( int nrl, int nrh, int ncl, int nch);
void    free_imatrix( int **m, int nrl, int nrh, int ncl, int nch);

float   *fvector( int nl, int nh);
void    free_fvector( float *v, int nl, int nh);
float   **fmatrix( int nrl, int nrh, int ncl, int nch);
void    free_fmatrix( float **m, int nrl, int nrh, int ncl, int nch);
float   **convert_fmatrix( float *m, int nrl, int nrh, int ncl, int nch);
void    free_convert_fmatrix( float **m, int nrl, int nrh, int ncl, int nch);

double  *dvector( int nl, int nh);
void    free_dvector( double *v, int nl, int nh);
double  **dmatrix( int nrl, int nrh, int ncl, int nch);
void    free_dmatrix( double **m, int nrl, int nrh, int ncl, int nch);
double  **convert_dmatrix( double *m, int nrl, int nrh, int ncl, int nch);
void    free_convert_dmatrix( double **m, int nrl, int nrh, int ncl, int nch);


/*
 *    dcomplex.c
 *
 *    Complex arithmetic and analysis in double precision.
 */
dcomplex *dcvector( int nl, int nh );
void  free_dcvector( dcomplex *v, int nl, int nh );
dcomplex **dcmatrix( int nrl, int nrh, int ncl, int nch );
void  free_dcmatrix( dcomplex **m, int nrl, int nrh, int ncl, int nch );

dcomplex add_dc( dcomplex a, dcomplex b );
dcomplex sub_dc( dcomplex a, dcomplex b );
dcomplex mul_dc( dcomplex a, dcomplex b );
dcomplex div_dc( dcomplex a, dcomplex b );
dcomplex mk_dc( double a, double b );
double abs_dc( dcomplex a );
double phase_dc( dcomplex a );
dcomplex conj_dc( dcomplex a );
dcomplex sqrt_dc( dcomplex a );
dcomplex iexp_dc( double a );
dcomplex exp_dc( dcomplex a );
dcomplex sin_dc( dcomplex a );
dcomplex cos_dc( dcomplex a );
dcomplex tan_dc( dcomplex a );
dcomplex sinh_dc( dcomplex a );
dcomplex cosh_dc( dcomplex a );
dcomplex tanh_dc( dcomplex a );
dcomplex log_dc( dcomplex a );
dcomplex log10_dc( dcomplex a );
dcomplex rmul_dc( double a, dcomplex b );
dcomplex imul_dc( double a, dcomplex b );


/*
 *    fcomplex.c
 *
 *    Complex arithmetic and analysis in single precision.
 */
fcomplex *fcvector( int nl, int nh );
void  free_fcvector( fcomplex *v, int nl, int nh );
fcomplex **fcmatrix( int nrl, int nrh, int ncl, int nch );
void  free_fcmatrix( fcomplex **m, int nrl, int nrh, int ncl, int nch );

fcomplex add_fc( fcomplex a, fcomplex b );
fcomplex sub_fc( fcomplex a, fcomplex b );
fcomplex mul_fc( fcomplex a, fcomplex b );
fcomplex div_fc( fcomplex a, fcomplex b );
fcomplex mk_fc( float a, float b );
float abs_fc( fcomplex a );
float phase_fc( fcomplex a );
fcomplex conj_fc( fcomplex a );
fcomplex sqrt_fc( fcomplex a );
fcomplex iexp_fc( float a );
fcomplex exp_fc( fcomplex a );
fcomplex sin_fc( fcomplex a );
fcomplex cos_fc( fcomplex a );
fcomplex tan_fc( fcomplex a );
fcomplex sinh_fc( fcomplex a );
fcomplex cosh_fc( fcomplex a );
fcomplex tanh_fc( fcomplex a );
fcomplex log_fc( fcomplex a );
fcomplex log10_fc( fcomplex a );
fcomplex rmul_fc( float a, fcomplex b );
fcomplex imul_fc( float a, fcomplex b );


/*
 *    dvector.c
 */
dmatrix32 cart_vec3vec2_dv( dvector3 U, dvector2 V );
dmatrix23 cart_vec2vec3_dv( dvector2 U, dvector3 V );
dvector3 mat32vec2_dv( dmatrix32 A, dvector2 V );
dvector2 mat23vec3_dv( dmatrix23 A, dvector3 V );
dmatrix32 trans_mat23_dv( dmatrix23 A );
dmatrix23 trans_mat32_dv( dmatrix32 A );
dmatrix3 mat32mat23_dv( dmatrix32 A, dmatrix23 B );
dmatrix32 mat33mat32_dv( dmatrix3 A, dmatrix32 B );
dmatrix32 mat32mat22_dv( dmatrix32 A, dmatrix2 B );
dmatrix2 mat23mat32_dv( dmatrix23 A, dmatrix32 B );
dmatrix23 mat22mat23_dv( dmatrix2 A, dmatrix23 B );
dmatrix23 mat23mat33_dv( dmatrix23 A, dmatrix3 B );

/*
 *    Algebra and analysis of double-precision three-dimensional vectors and matrices.
 */

dvector3 *dv3vector( int nlo, int nhi );
void  free_dv3vector( dvector3 *v, int nl, int nh );
dvector2 *dv2vector( int nlo, int nhi );
void  free_dv2vector( dvector2 *v, int nl, int nh );

double *convert_dvector_dv3( dvector3 a );
void free_convert_dvector_dv3( double *v );
double **convert_dmatrix_dv3( dmatrix3 a );
void free_convert_dmatrix_dv3( double **m );
dvector3 mk_dvector_dv3( double *v );
dmatrix3 mk_dmatrix_dv3( double **m );
dvector3 mk_dv3( double x, double y, double z );
dvector3 add_dv3( dvector3 a, dvector3 b );
dvector3 sub_dv3( dvector3 a, dvector3 b );
double mod_dv3( dvector3 a );
dvector3 norm_dv3( dvector3 a );
double dot_dv3( dvector3 a, dvector3 b );
double az_dv3( dvector3 a );
double el_dv3( dvector3 a );
dvector3 cross_dv3( dvector3 a, dvector3 b );
dvector3 rmul_dv3( double a, dvector3 b );
dmatrix3 cart_dv3( dvector3 a, dvector3 b );
double trace_dv3( dmatrix3 A );
double det_dv3( dmatrix3 A );
dmatrix3 inv_dv3( dmatrix3 A );
dmatrix3 matmat_dv3( dmatrix3 A, dmatrix3 B );
dvector3 matvec_dv3( dmatrix3 A, dvector3 X );
dvector3 vecmat_dv3( dvector3 X, dmatrix3 A );
double vecmatvec_dv3( dvector3 Y, dmatrix3 A, dvector3 X);
dmatrix3 trans_dv3( dmatrix3 A );
/*
 *    Algebra and analysis of double-precision two-dimensional vectors and matrices.
 */
double *convert_dvector_dv2( dvector2 a );
void free_convert_dvector_dv2( double *v );
double **convert_dmatrix_dv2( dmatrix2 a );
void free_convert_dmatrix_dv2( double **m );
dvector2 mk_dvector_dv2( double *v );
dmatrix2 mk_dmatrix_dv2( double **m );
dvector2 mk_dv2( double x, double y );
dvector2 add_dv2( dvector2 a, dvector2 b );
dvector2 sub_dv2( dvector2 a, dvector2 b );
double mod_dv2( dvector2 a );
dvector2 norm_dv2( dvector2 a );
double dot_dv2( dvector2 a, dvector2 b );
double az_dv2( dvector2 a );
dvector2 rmul_dv2( double a, dvector2 b );
dmatrix2 cart_dv2( dvector2 a, dvector2 b );
double trace_dv2( dmatrix2 A );
double det_dv2( dmatrix2 A );
dmatrix2 inv_dv2( dmatrix2 A );
dmatrix2 matmat_dv2( dmatrix2 A, dmatrix2 B );
dvector2 matvec_dv2( dmatrix2 A, dvector2 X );
dvector2 vecmat_dv2( dvector2 X, dmatrix2 A );
double vecmatvec_dv2( dvector2 Y, dmatrix2 A, dvector2 X);
dmatrix2 trans_dv2( dmatrix2 A );



/*
 *    fvector.c
 */
fmatrix32 cart_vec3vec2_fv( fvector3 U, fvector2 V );
fmatrix23 cart_vec2vec3_fv( fvector2 U, fvector3 V );
fvector3 mat32vec2_fv( fmatrix32 A, fvector2 V );
fvector2 mat23vec3_fv( fmatrix23 A, fvector3 V );
fmatrix32 trans_mat23_fv( fmatrix23 A );
fmatrix23 trans_mat32_fv( fmatrix32 A );
fmatrix3 mat32mat23_fv( fmatrix32 A, fmatrix23 B );
fmatrix32 mat33mat32_fv( fmatrix3 A, fmatrix32 B );
fmatrix32 mat32mat22_fv( fmatrix32 A, fmatrix2 B );
fmatrix2 mat23mat32_fv( fmatrix23 A, fmatrix32 B );
fmatrix23 mat22mat23_fv( fmatrix2 A, fmatrix23 B );
fmatrix23 mat23mat33_fv( fmatrix23 A, fmatrix3 B );
/*
 *    Algebra and analysis of single-precision three-dimensional vectors and tensors.
 */
fvector3 *fv3vector( int nlo, int nhi );
void  free_fv3vector( fvector3 *v, int nl, int nh );
fvector2 *fv2vector( int nlo, int nhi );
void  free_fv2vector( fvector2 *v, int nl, int nh );

float *convert_fvector_fv3( fvector3 a );
void free_convert_fvector_fv3( float *v );
float **convert_fmatrix_fv3( fmatrix3 a );
void free_convert_fmatrix_fv3( float **m );
fvector3 mk_fvector_fv3( float *v );
fmatrix3 mk_fmatrix_fv3( float **m );
fvector3 mk_fv3( float x, float y, float z );
fvector3 add_fv3( fvector3 a, fvector3 b );
fvector3 sub_fv3( fvector3 a, fvector3 b );
float mod_fv3( fvector3 a );
fvector3 norm_fv3( fvector3 a );
float dot_fv3( fvector3 a, fvector3 b );
float az_fv3( fvector3 a );
float el_fv3( fvector3 a );
fvector3 cross_fv3( fvector3 a, fvector3 b );
fvector3 rmul_fv3( float a, fvector3 b );
fmatrix3 cart_fv3( fvector3 a, fvector3 b );
float trace_fv3( fmatrix3 A );
float det_fv3( fmatrix3 A );
fmatrix3 inv_fv3( fmatrix3 A );
fmatrix3 matmat_fv3( fmatrix3 A, fmatrix3 B );
fvector3 matvec_fv3( fmatrix3 A, fvector3 X );
fvector3 vecmat_fv3( fvector3 X, fmatrix3 A );
float vecmatvec_fv3( fvector3 Y, fmatrix3 A, fvector3 X);
fmatrix3 trans_fv3( fmatrix3 A );
/*
 *    Algebra and analysis of single-precision two-dimensional vectors and tensors.
 */
float *convert_fvector_fv2( fvector2 a );
void free_convert_fvector_fv2( float *v );
float **convert_fmatrix_fv2( fmatrix2 a );
void free_convert_fmatrix_fv2( float **m );
fvector2 mk_fvector_fv2( float *v );
fmatrix2 mk_fmatrix_fv2( float **m );
fvector2 mk_fv2( float x, float y );
fvector2 add_fv2( fvector2 a, fvector2 b );
fvector2 sub_fv2( fvector2 a, fvector2 b );
float mod_fv2( fvector2 a );
fvector2 norm_fv2( fvector2 a );
float dot_fv2( fvector2 a, fvector2 b );
float az_fv2( fvector2 a );
fvector2 rmul_fv2( float a, fvector2 b );
fmatrix2 cart_fv2( fvector2 a, fvector2 b );
float trace_fv2( fmatrix2 A );
float det_fv2( fmatrix2 A );
fmatrix2 inv_fv2( fmatrix2 A );
fmatrix2 matmat_fv2( fmatrix2 A, fmatrix2 B );
fvector2 matvec_fv2( fmatrix2 A, fvector2 X );
fvector2 vecmat_fv2( fvector2 X, fmatrix2 A );
float vecmatvec_fv2( fvector2 Y, fmatrix2 A, fvector2 X);
fmatrix2 trans_fv2( fmatrix2 A );


/*
 *    dlinear.c
 *
 *    Linear algebra in double precision for arbitrary dimension matrices and vectors.
 */
void dcopyvec( double *a, int nrowlo, int nrowhi, double *b );
void dcopymat( double **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **b );
void dmodvec( double *a, int nlo, int nhi, double *c );
void dcartvecvec( double *a, double *b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c );
void dvecvec( double *a, double *b, int nlo, int nhi, double *c );
void dmatvec( double **a, double *b, int nrowlo, int nrowhi, int nlo, int nhi, double *c );
void dvecmat( double *b, double **a, int nlo, int nhi, int ncollo, int ncolhi, double *c );
void dvecmatvec( double *v, double **a, double *u, int n1lo, int n1hi, int n2lo, int n2hi, double *c );
void dmatmat( double **a, double **b, int nrowlo, int nrowhi, int nlo, int nhi, int ncollo, int ncolhi, double **c );
void daddvec( double *a, double *b, int nrowlo, int nrowhi, double *c );
void dsubvec( double *a, double *b, int nrowlo, int nrowhi, double *c );
void dmulvec( double a, double *b, int nrowlo, int nrowhi, double *c );
void daddmat( double **a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c );
void dsubmat( double **a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c );
void dmulmat( double a, double **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **c );
void dtransmat( double **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, double **atrans );
void dsymmat( double **a, int nlo, int nhi, double **asymm, double **aanti );
void dtracemat( double **a, int nlo, int nhi, double *trace );
void dgramm_schmidt( double **a, int nlo, int nhi, double **aortho );


/*
 *    flinear.c
 *
 *    Linear algebra in single precision for arbitrary dimension matrices and vectors.
 */
void fcopyvec( float *a, int nrowlo, int nrowhi, float *b );
void fcopymat( float **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **b );
void fmodvec( float *a, int nlo, int nhi, float *c );
void fcartvecvec( float *a, float *b, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **c );
void fvecvec( float *a, float *b, int nlo, int nhi, float *c );
void fmatvec( float **a, float *b, int nrowlo, int nrowhi, int nlo, int nhi, float *c );
void fvecmat( float *b, float **a, int nlo, int nhi, int ncollo, int ncolhi, float *c );
void fvecmatvec( float *v, float **a, float *u, int n1lo, int n1hi, int n2lo, int n2hi, float *c );
void fmatmat( float **a, float **b, int nrowlo, int nrowhi, int nlo, int nhi, int ncollo, int ncolhi, float **c );
void faddvec( float *a, float *b, int nrowlo, int nrowhi, float *c );
void fsubvec( float *a, float *b, int nrowlo, int nrowhi, float *c );
void fmulvec( float a, float *b, int nrowlo, int nrowhi, float *c );
void faddmat( float **a, float **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **c );
void fsubmat( float **a, float **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **c );
void fmulmat( float a, float **b, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **c );
void ftransmat( float **a, int nrowlo, int nrowhi, int ncollo, int ncolhi, float **atrans );
void fsymmat( float **a, int nlo, int nhi, float **asymm, float **aanti );
void ftracemat( float **a, int nlo, int nhi, float *trace );
void fgramm_schmidt( float **a, int nlo, int nhi, float **aortho );


/*
 *    roots.c
 */
int root_bracket( double (*funk)(double), double *x1, double *x2 );
void root_brackets( double (*funk)(double),
                    double x1, double x2, int n,
                    double *xb1, double *xb2, int *nb );
double root_bisect( double (*funk)(double), double x1, double x2, double err );
double root_secant( double (*funk)(double), double x1, double x2, double err );
double root_false_position( double (*funk)(double), double x1, double x2, double err );
double root_ridder( double (*funk)(double), double x1, double x2, double err );
double root_brent( double (*funk)(double), double x1, double x2, double err );
double root_newton( void (*funky)( double, double *, double * ), double x1, double x2, double err );
double root_safe( void (*funky)( double, double *, double * ), double x1, double x2, double err );


double min_golden( double ax, double bx, double cx, double (*funk)(double), double tol, double *xmin );
double min_brent( double ax, double bx, double cx, double (*funk)(double), double tol, double *xmin );


/*
 *    sorting.c
 */
void sort_pick( int n, double *arr );
void sort_pick2( int n, double *arr, double *brr );
void sort_shell( unsigned long n, double *a );


/*
 *    interpolation.c
 */
void poly_interp( double *xa, double *ya, int n, double x, double *y, double *dy );
void rat_interp( double *xa, double *ya, int n, double x, double *y, double *dy );
void spline( double *x, double *y, int n, double yp_1, double yp_n, double *ypp );
double spline_interp( double *xa, double *ya, double *yapp, int na, double x );
double linear_interp( double *xa, double *ya, int na, double x );


/*
 *    chebyshev.c
 */
void chebyshev_fit( double a, double b, double *c, int n, double (*funk)( double ) );
double chebyshev_eval( double a, double b, double *c, int m, double x );
void chebyshev_deriv( double a, double b, double *c, double *cder, int m );
void chebyshev_integ( double a, double b, double *c, double *cint, int m);
void chebyshev_2_poly( double a, double b, double *c, double *d, int m );
void poly_2_chebyshev( double a, double b, double *d, double *c, int m );
void poly_economize( double a, double b, double *d, double *c, int m, int *mfew, double err );
double clenshaw_curtis_quad( double a, double b, double *c, int m, double err );


/*
 *    poly.c
 */
void poly_shift_coeff( double a, double b, double *c, int n );


/*
 *    gauss_jordan.c
 */
void gauss_jordan( double **a, int n, double **b, int m );


/*
 *    lu_decomp.c
 */
void lu_decomp( double **a, int n, int *index, double *parity);
void lu_backsub( double **a, int n, int *index, double *b);
void lu_improve( double **a, double **alud, int n, int *index, double *b, double *x);
void lu_invert( double **alud, int n, int *index, double **ainv );
double lu_determinant( double **alud, int n, double parity );
double lu_trace( double **alud, int n );


/*
 *    tridiag.c
 */
void tridiagonal( double *a, double *b, double *c, double *r, double *u, int n);
void cyclic( double *a, double *b, double *c, double alpha, double beta, double *r, double *x, unsigned long n );


/*
 *    sv_decomp.c
 */
void sv_decomp( double **aa, int m, int n, double *w, double **v );
void sv_backsub( double **u, double *w, double **v, int m, int n, double *b, double *x );
void sv_improve( double **a, double **u, double *w, double **v, int m, int n, double *b, double *x );


/*
 *    cholesky_decomp.c
 */
void cholesky_decomp( double **a, int n, double *d );
void cholesky_backsub( double **a, int n, double *d, double *b, double *x );
void cholesky_invert( double **a, int n, double *d );


/*
 *    qr_decomp.c
 */
void qr_decomp( double **a, int n, double *c, double *d, int *sing );
void qr_backsub( double **a, int n, double *c, double *d, double *b );
void r_backsub( double **a, int n, double *d, double *b );
void qr_update( double **r, double **qt, int n, double *u, double *v );
void jacobi_rotate( double **r, double **qt, int n, int i, double a, double b );


/*
 *    fit.c
 */
void fit( double *x, double *y, int ndata, double *sigma, int mwt,
          double *a, double *b, double *sigmaa,
          double *sigmab, double *chisqr, double *q );
void svd_fit( double *x, double *y, double *sig, int ndata, double *a, int ma,
              double **u, double **v, double *w, double *tol, double *chisqr,
              void (*funcs)( double, double *, int) );
void svd_var( double **v, int ma, double *w, double **cov );
void funcs_poly( double x, double *afunc, int n );
void funcs_legendre( double x, double *afunc, int n );



/*
 *    gamma.c
 */
double ln_gamma( double x);
double factorial( int n);
double bin_coeff( int n, int k);
double ln_factorial( int n);
double beta( double z, double w);
double gamma_p( double a, double x);
double gamma_q( double a, double x);
void   gamma_series( double *gamser, double a, double x, double *lngam);
void   gamma_cont_fraction( double *gamconfrac,
                             double a, double x, double *lngam);
double error_func( double x);
double error_func_comp( double x);


/*
 *    sph_harm.c
 */
double spherical_harmonic( int l, int m, double theta, double phi);
double legendre_poly( int l, int m, double x);


/*
 *    poly_jacobi.c
 */
double jacobi_poly( int k, double alpha, double beta, double x);


/*
 *    poly_gegenbauer.c
 */
double gegenbauer_poly( int n, double alpha, double x);


/*
 *    poly_laguerre.c
 */
double laguerre_poly( int n, double alpha, double x);


/*
 *    poly_hermite.c
 */
double hermite_h0( double x);
double hermite_h1( double x);
double hermite_h( int n, double x);


/*
 *    poly_legendre.c
 */
double legendre_p0( double x);
double legendre_p1( double x);
double legendre_p( int l, double x);
double legendre_q0( double x);
double legendre_q1( double x);
double legendre_q( int l, double x);


/*
 *    poly_chebyshev.c
 */
double chebyshev_t0( double x );
double chebyshev_t1( double x );
double chebyshev_t( int n, double x );


/*
 *    bessel_general.c
 */
void bessel_jy( double xnu, double x, double *rj, double *ry, double *rjp, double *ryp );
void bessel_ik( double xnu, double x, double *ri, double *rk, double *rip, double *rkp );
void bessel_cheb( double x, double *gam1, double *gam2, double *gampl, double *gammi );
void airy( double x, double *ai, double *bi, double *aip, double *bip );
void wairy( double x, dcomplex *w1, dcomplex *w2, dcomplex *w1p, dcomplex *w2p );
void airy_root( double *alpha, double *beta, double *alpha_p, double *beta_p, int s );
void sph_bessel( int n, double x, double *sj, double *sy, double *sjp, double *syp );


/*
 *    bessel_integer.c
 */
double bessel_j0( double x );
double bessel_j1( double x );
double bessel_j( int n, double x );
double bessel_y0( double x );
double bessel_y1( double x );
double bessel_y( int n, double x );
double bessel_i0( double x );
double bessel_i1( double x );
double bessel_i( int n, double x );
double bessel_k0( double x );
double bessel_k1( double x );
double bessel_k( int n, double x );


/*
 *    bessel_spherical.c
 */
double sph_bessel_j0( double x);
double sph_bessel_j1( double x);
double sph_bessel_j( int n, double x);

double sph_bessel_y0( double x);
double sph_bessel_y1( double x);
double sph_bessel_y( int n, double x);

double sph_bessel_i0( double x);
double sph_bessel_i1( double x);
double sph_bessel_i( int n, double x);

double sph_bessel_k0( double x);
double sph_bessel_k1( double x);
double sph_bessel_k( int n, double x);


/*
 *    int_elliptic.c
 */
double carlson_rf( double x, double y, double z );
double carlson_rd( double x, double y, double z );
double carlson_rj( double x, double y, double z, double p );
double carlson_rc( double x, double y );
double legendre_e( double phi, double ak );
double legendre_f( double phi, double ak );
double legendre_pi( double phi, double en, double ak );
void jacobian_sncndn( double uu, double emmc, double *sn, double *cn, double *dn );


/*
 *    int_exp.c
 */
double exp_int( int n, double x );
double ei( double x );


/*
 *    int_cossin.c
 */
void cisi( double x, double *ci, double *si );


/*
 *    int_dawson.c
 */
double dawson( double x );


/*
 *    int_fresnel.c
 */
double fresnel_c( double x );
double fresnel_s( double x );
double fresnel( double x, double *c, double *s );
dcomplex fock( double x );
double fock_r( double x );
double fock_i( double x );



/*
 *    integration.c
 */
double trapezoid( double (*funk)(double), double a, double b, int n );
double midpoint( double (*funk)(double), double a, double b, int n );
double midpoint_inv( double (*funk)(double), double a, double b, int n );
double midpoint_inv_sqrt_lower( double (*funk)(double), double a, double b, int n );
double midpoint_inv_sqrt_upper( double (*funk)(double), double a, double b, int n );
double midpoint_exp( double (*funk)(double), double a, double b, int n );
double dumb_trapezoid( double (*funk)(double), double a, double b, int n );
double quad_trapezoid( double (*funk)(double), double a, double b, double err );
double dumb_simpson( double (*funk)(double), double a, double b, int n );
double quad_simpson( double (*funk)(double), double a, double b, double err );
double dumb_romberg( double (*funk)(double), double a, double b, int n );
double quad_romberg( double (*funk)(double), double a, double b, double err );
double quad_romberg_open( double (*funk)(double), double a, double b, double err,
                       double (*choose)( double (*)(double), double, double, int ) );



/*
 *    diff.c
 */
double diff_ridder( double (*funk)(double), double x, double h, double *err );



/*
 *    gauss_quad.c
 */
double quad_gauss_legendre( double (*funk)(double), double *x, double *w, int n );
double quad_gauss( double (*funk)(double), double a, double b, double *x, double *w, int n );
void gauss_legendre( int n, double xa, double xb, double *x, double *w );
void gauss_laguerre( double *x, double *w, int n, double alpha );
void gauss_hermite( double *x, double *w, int n );
void gauss_jacobi( double *x, double *w, int n, double alpha, double beta );
void gauss_chebyshev( double *x, double *w, int n );
double gauss_crap( double (*funk)(double), double a, double b, int n );
double dumb_gauss_crap( double (*funk)(double), double a, double b, int n );
double quad_gauss_crap( double (*funk)(double), double a, double b, double err );



/*
 *    ode.c
 */
void runge_kutta_4( double *y, double *dydx, int n,
                    double x, double h, double *yout,
                    void (*derivs)( double, double *, double *) );

void dumb_runge_kutta( double *ystart, int nvar, double x1, double x2, int nstep,
                       double *runge_kutta_xx, double **runge_kutta_yy,
                       void (*derivs)( double, double *, double *) );

void quad_runge_kutta( double *y, double *dydx, int n, double *x, double htry,
                       double err, double *yscale, double *hdid, double *hnext,
                       void (*derivs)( double, double *, double *) );

void cash_karp_rk( double *y, double *dydx, int n,
                   double x, double h, double *yout, double *yerr,
                   void (*derivs)( double, double *, double *) );

void quad_cash_karp_rk( double *y, double *dydx, int n, double *x,
                        double htry, double err, double *yscale,
                        double *hdid, double *hnext,
                        void (*derivs)( double, double *, double *) );

void modified_midpoint( double *y, double *dydx, int nvar, double xs,
                        double htot, int nstep, double *yout,
                        void (*derivs)( double, double *, double *) );

void stoermer( double *y, double *d2y, int nv, double xs,
               double htot, int nstep, double *yout,
               void (*derivs)( double, double *, double *) );

void bulirsch_stoer( double *y, double *dydx, int nv, double *xx,
                     double htry, double err, double *yscale,
                     double *hdid, double *hnext,
                     void (*derivs)( double, double *, double *) );

void ode_integrate( double *ystart, int nvar, double x1, double x2,
                    double err, double h1, double hmin, 
                    int *nok, int *nbad, 
                    void (*derivs)( double, double *, double * ), 
                    void (*stepper)( double *, double *, int, double *, double ,
                                     double, double *, double *, double *,
                                     void (*)( double, double *, double *) ) );

/*
 *    fourier.c
 */
void fourier1( double *data, unsigned long nn, int isign );
void four_real( double *data, unsigned long nn, int isign );
void four_sin( double *y, int n );
void four_cos1( double *y, int n );
void four_cos2( double *y, int n );



/*
 *    quadratic.c
 */
int quadratic( double a0, double a1, double a2, double *r1, double *r2 );
void quad_complex( dcomplex a0, dcomplex a1, dcomplex a2, dcomplex *r1, dcomplex *r2 );


/*
 *    util.c
 */

double dsgn( double a );
float fsgn( float a );
int isgn( int a );

double dsign( double a, double b );
float fsign( float a, float b );
int isign( int a, int b );

double dmin( double a, double b );
//float fmin( float a, float b );
int imin( int a, int b );

double dmin3( double a, double b, double c );
float fmin3( float a, float b, float c );
int imin3( int a, int b, int c );

double dmax( double a, double b );
//float fmax( float a, float b );
int imax( int a, int b );

double dmax3( double a, double b, double c );
float fmax3( float a, float b, float c );
int imax3( int a, int b, int c );

void dswap( double *a, double *b );
void fswap( float *a, float *b );
void iswap( int *a, int *b );

void drpermute3( double *a, double *b, double *c );
void frpermute3( float *a, float *b, float *c );
void irpermute3( int *a, int *b, int *c );

void dlpermute3( double *a, double *b, double *c );
void flpermute3( float *a, float *b, float *c );
void ilpermute3( int *a, int *b, int *c );

void drshift3( double *a, double *b, double *c );
void frshift3( float *a, float *b, float *c );
void irshift3( int *a, int *b, int *c );

void dlshift3( double *a, double *b, double *c );
void flshift3( float *a, float *b, float *c );
void ilshift3( int *a, int *b, int *c );

void drshift4( double *a, double *b, double *c, double *d);
void frshift4( float *a, float *b, float *c, float *d);
void irshift4( int *a, int *b, int *c, int *d);

void dlshift4( double *a, double *b, double *c, double *d);
void flshift4( float *a, float *b, float *c, float *d);
void ilshift4( int *a, int *b, int *c, int *d);

double dsqr( double a );
float fsqr( float a );
int isqr( int a );

double dcub( double a );
float fcub( float a );
int icub( int a );

double dpow4( double a );
float fpow4( float a );
int ipow4( int a );

double dpow5( double a );
float fpow5( float a );
int ipow5( int a );

double dpythag( double a, double b );
float fpythag( float a, float b );



/*
 *    wave_function.c
 */
double  radial_Hydrogen( int Z, int n, int l, double r, double *E_n);
double  harmonic( double mu, double omega, int n, double x, double *E_n);
double  radial_harmonic_2d( double mu, double omega, int n, int m, double p, double *E_nm);
double  radial_harmonic_3d( double mu, double omega, int n, int l, double r, double *E_nl);


#ifdef __cplusplus
}
#endif

#endif  /*  __NRIC_H__  */

