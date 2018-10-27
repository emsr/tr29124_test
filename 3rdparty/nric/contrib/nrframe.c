/* nrframe.c **************************************************************

  NR Library Package - Homogeneous Transformation Frames

  James Trevelyan,  University of Western Australia
  Revision 2   January 1996

**************************************************************************/

/* define V_CHECK to include validity checks on matrices */
#define V_CHECK

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrlin.h"

#define True 1
#define False 0
#define X 1
#define Y 2
#define Z 3
#define U 4

#ifdef TEST
extern FILE *outfile;
extern double **F1;
#endif

double **SCR;

double **dframe( void )
{
   return ( dmatrix( 1, 4, 1, 4 ) );
}

double *d4vector( void )
{
   return ( dvector( 1, 4 ) );
}

void chainmult( double **a, double **b )
/* multiply 4*4 matrices a, b into a, using SCR scratch matrix above */
{   
   static int once = 0;

   if ( !once ) {
      SCR = dmatrix( 1,4, 1,4 );
      once = True;
   }
#ifdef V_CHECK
   if ( !valid_dframe( a ) ) 
      nrerror("Invalid 1st matrix: chainmult\n"); 
   if ( !valid_dframe( b ) ) 
      nrerror("Invalid 2nd matrix: chainmult\n"); 
#endif
   dmmult( a, 4, 4, b, 4, 4, SCR );
   dmcopy( SCR, 4, 4, a );
}


void rotframe( double **a, int axis, double theta) /* make a rotation transform */
{
   double ca, sa;
   int i,j;

#ifdef V_CHECK
   if ( !valid_dframe( a ) ) 
      nrerror("Invalid matrix: rotframe\n");
#endif
   ca = cos(theta);
   sa = sin(theta);
   for (i=1; i<=4; i++)
      for (j=1; j<=4; j++) a[i][j] = 0.0;
   switch( axis) {
     case 3:
      a[4][4] = a[3][3] = 1.0;
      a[2][2] = a[1][1] = ca;
      a[1][2] = -sa;
      a[2][1] = sa;
      break;
     case 1:
      a[4][4] = a[1][1] = 1.0;
      a[2][2] = a[3][3] = ca;
      a[2][3] = -sa;
      a[3][2] = sa;
      break;
     case 2:
      a[4][4] = a[2][2] = 1.0;
      a[1][1] = a[3][3] = ca;
      a[3][1] = -sa;
      a[1][3] = sa;
      break;
     default: 
      fprintf(stderr, "rotframe: axis not X, Y or Z\n");
      exit(1);
      break;
   }
}

void transframe( double **a, double dx, double dy, double dz) 
/* make a translation tr'm */
{
   int i,j;

#ifdef V_CHECK
   if ( !valid_dframe( a ) ) 
      nrerror("Invalid matrix: transframe\n");    
#endif
   for (i=1; i<=4; i++)
      for (j=1; j<=4; j++) a[i][j] = 0.0;
   a[4][4] = a[3][3] = 1.0;
   a[2][2] = a[1][1] = 1.0;
   a[1][4] = dx;
   a[2][4] = dy;
   a[3][4] = dz;
}


void orthog( double **a ) /* orthogonalize a frame */
{
   int i, j, k;
   double s;

#ifdef V_CHECK
   if ( !valid_dframe( a ) ) 
      nrerror("Invalid matrix: orthog\n");    
#endif
   /* cross j, k to produce i */
   a[1][1] = a[2][2]*a[3][3] - a[3][2]*a[2][3];
   a[2][1] = a[3][2]*a[1][3] - a[1][2]*a[3][3];
   a[3][1] = a[1][2]*a[2][3] - a[2][2]*a[1][3];

   /* normalize i, k */
   s = sqrt(DSQR(a[1][3]) + DSQR(a[2][3]) + DSQR(a[3][3]));
   for (i=1; i<=3; i++) a[i][3] /= s;
   s = sqrt(DSQR(a[1][1]) + DSQR(a[2][1]) + DSQR(a[3][1]));
   for (i=1; i<=3; i++) a[i][1] /= s;

   /* cross k, i to produce j */
   a[1][2] = a[2][3]*a[3][1] - a[3][3]*a[2][1];
   a[2][2] = a[3][3]*a[1][1] - a[1][3]*a[3][1];
   a[3][2] = a[1][3]*a[2][1] - a[2][3]*a[1][1];

}

void orthogKI( double **a ) /* orthogonalize a frame using K & I vector */
{
   int i, j, k;
   double s;

#ifdef V_CHECK
   if ( !valid_dframe( a ) )
      nrerror("Invalid matrix: orthog\n");
#endif
   /* cross k, i to produce j */
   a[1][2] = a[2][3]*a[3][1] - a[3][3]*a[2][1];
   a[2][2] = a[3][3]*a[1][1] - a[1][3]*a[3][1];
   a[3][2] = a[1][3]*a[2][1] - a[2][3]*a[1][1];

   /* normalize j, k */
   s = sqrt(DSQR(a[1][3]) + DSQR(a[2][3]) + DSQR(a[3][3]));
   for (i=1; i<=3; i++) a[i][3] /= s;
   s = sqrt(DSQR(a[1][2]) + DSQR(a[2][2]) + DSQR(a[3][2]));
   for (i=1; i<=3; i++) a[i][2] /= s;

   /* cross j, k to produce i */
   a[1][1] = a[2][2]*a[3][3] - a[3][2]*a[2][3];
   a[2][1] = a[3][2]*a[1][3] - a[1][2]*a[3][3];
   a[3][1] = a[1][2]*a[2][3] - a[2][2]*a[1][3];

}

void dframeinverse( double **E, double **Ein ) /* form inverse frame */
{
   double a[6], b[6];
   int i, j;

   validate_dvector( a, 1, 4);
   validate_dvector( b, 1, 4);

   for (i=1; i<=3; i++) {
      for (j=1; j<=3; j++) {
         Ein[i][j] = E[j][i];
      }
      Ein[4][i] = 0.0;
      Ein[i][4] = 0.0;
   }
   Ein[4][4] = 1.0;

   dgetcolumn(E, 4, a, 4);
   dmvmult(Ein, 4, 4, a, 4, b);
   for (i=1; i<=3; i++) Ein[i][4] = -b[i];
}


void dcross4 ( double *a, double *b, double *c ) /* cross two 4-vectors */
{
   double x, y, z;

#ifdef V_CHECK
   if ( !valid_d4vector_b( a ) ) 
      nrerror("Invalid vector 1: dcross4\n"); 
   if ( !valid_d4vector_b( b ) ) 
      nrerror("Invalid vector 2: dcross4\n"); 
   if ( !valid_d4vector_b( c ) ) 
      nrerror("Invalid vector 3: dcross4\n"); 
#endif
   x = a[2]*b[3] - a[3]*b[2];
   y = a[3]*b[1] - a[1]*b[3];
   z = a[1]*b[2] - a[2]*b[1];
   c[1] = x; 
   c[2] = y; 
   c[3] = z; 
   c[4] = 0.0;
}

double column_dot( double **G1, int col1, double **G2, int col2 )
/* form vector dot product of two columns (1st 3 elements only) */
{
   int i;
   double sum;

#ifdef V_CHECK
   if ( !valid_dframe( G1 ) ) 
      nrerror("Invalid 1st frame: column_dot\n"); 
   if ( !valid_dframe( G2 ) ) 
      nrerror("Invalid 2nd frame: column_dot\n"); 
#endif
   sum = 0.0;
   for ( i=1; i<=3; i++ ) sum += G1[i][col1]*G2[i][col2];
   return(sum);
}

void make_d_omega( double **G1, double **G2, double *w)
/* find small rotational difference between two frames G1, G2 */
/* ref: Shear Magic 4A.27 p129 */
{
#ifdef V_CHECK
   if ( !valid_dframe( G1 ) ) 
      nrerror("Invalid 1st frame: make_d_omega\n");
   if ( !valid_dframe( G2 ) ) 
      nrerror("Invalid 2nd frame: make_d_omega\n");   
#endif
   w[1] = column_dot( G1, 3, G2, 2 );
   w[2] = column_dot( G1, 1, G2, 3 );
   w[3] = column_dot( G1, 2, G2, 1 );
   w[4] = 0.0;
}

void compare_link_frames( double **G1, double **G2, double *wp, double *wr, double *ww )
{
   double wrot[6], z1[6], z2[6], zz[6];

   validate_dvector(wrot, 1, 4);
   validate_dvector(z1, 1, 4);
   validate_dvector(z2, 1, 4);
   validate_dvector(zz, 1, 4);

#ifdef V_CHECK
   if ( !valid_dframe( G1 ) ) 
      nrerror("Invalid 1st frame: compare_link_frames\n");    
   if ( !valid_dframe( G2 ) ) 
      nrerror("Invalid 2nd frame: compare_link_frames\n");    
#endif
   *wp = *wr = *ww = 0.0;
   dgetcolumn( G1, Z, z1, 4 ); /* z axis comparison by cross product */
   dgetcolumn( G2, Z, z2, 4 );
   dcross4( z1, z2, zz );
   *wr = dvmag( zz, 3 );

   dgetcolumn( G1, U, z1, 4 ); /* origin comparison */
   dgetcolumn( G2, U, z2, 4 );
   dvsub( z1, 4, z2, zz );
   *wp = dvmag( zz, 3 );

   /*  printf("wp, wr, ww: %8.4lf %8.4lf %8.4lf\n", *wp, *wr, *ww ); */
}

void compare_frames( double **G1, double **G2, double *wp, double *wr, double *ww )
{
   double wrot[6], z1[6], z2[6], zz[6];

   validate_dvector(wrot, 1, 4);
   validate_dvector(z1, 1, 4);
   validate_dvector(z2, 1, 4);
   validate_dvector(zz, 1, 4);

#ifdef V_CHECK
   if ( !valid_dframe( G1 ) ) 
      nrerror("Invalid 1st frame: compare_frames\n"); 
   if ( !valid_dframe( G2 ) ) 
      nrerror("Invalid 2nd frame: compare_frames\n"); 
#endif
   *wp = *wr = *ww = 0.0;
   make_d_omega( G1, G2, wrot );   
   *wr = dvmag( wrot, 3 );
   *wp = sqrt( DSQR(G1[X][U] - G2[X][U]) 
      + DSQR(G1[Y][U] - G2[Y][U]) 
      + DSQR(G1[Z][U] - G2[Z][U]) );
   *ww = sqrt( DSQR(*wr) + DSQR(*wp) );
   /*  printf("wp, wr, ww: %8.4lf %8.4lf %8.4lf\n", *wp, *wr, *ww ); */
}

/*------------------------------------------------------------*/
void FindAxis(double **T1, double **T2, double *axis,
                double *sphi, double *cphi,int *twitching)
{
   /*  Routine to find the axis of rotation between two frames and the angle between them.
   */

   static double **FF;
   static int done_once = False;
   double ca[6],cb[6],cc[6],cd[6]; /*   workspace vectors   */
   double da[6],db[6],dc[6],dd[6]; /*   workspace vectors   */
   double d1, d2, d3, dmax;

   validate_dvector(ca, 1, 4);
   validate_dvector(cb, 1, 4);
   validate_dvector(cc, 1, 4);
   validate_dvector(cd, 1, 4);
   validate_dvector(da, 1, 4);
   validate_dvector(db, 1, 4);
   validate_dvector(dc, 1, 4);
   validate_dvector(dd, 1, 4);

   if (!done_once) {
      FF = dmatrix( 1, 4, 1, 4 );
      done_once = True;
   }

#ifdef V_CHECK
   if ( !valid_dframe( T1 ) ) 
      nrerror("Invalid 1st frame: FindAxis\n");   
   if ( !valid_dframe( T2 ) ) 
      nrerror("Invalid 2nd frame: FindAxis\n");   
   if ( !valid_d4vector_b( axis ) )
      nrerror("Invalid axis vector: FindAxis\n"); 
#endif

   *twitching = False;

   /*  Find axis of rotation */

   dmsub(T2, 4, 4, T1, FF);             /*  Difference of transforms    */
   dgetcolumn(FF, X, ca, 4);
   dgetcolumn(FF, Y, cb, 4);
   dgetcolumn(FF, Z, cc, 4);
   dcross4( ca, cb, da );
   dcross4( cb, cc, db );
   dcross4( cc, ca, dc );
   d1 = dvmnorm( da, 3 );
   d2 = dvmnorm( db, 3 );
   d3 = dvmnorm( dc, 3 );
   if ( d1 < d2 ) {
      if ( d2 < d3 ) {
         dvcopy( dc, 4, axis );
         dmax = d3;
      } 
      else {
         dvcopy( db, 4, axis );
         dmax = d2;
      }
   } 
   else {
      if ( d1 < d3 ) {
         dvcopy( dc, 4, axis );
         dmax = d3;
      } 
      else {
         dvcopy( da, 4, axis );
         dmax = d1;
      }
   }
   /*    printf("dmax = %12.6lf\n", dmax); */
   if( dmax < 0.000001 ) {                  /* there is no rotation, or we have a problem  */
      d1 = column_dot( T1, X, T2, X);
      d2 = column_dot( T1, Y, T2, Y);
      d3 = column_dot( T1, Z, T2, Z);
      if( d1<0.5 || d2<0.5 || d3<0.5 ) { /* problem */
         fprintf( stderr, "FindAxis: Frames opposed to each other\n");
         dgetcolumn(T1, Z, axis, 4);      /* Set axis to k vector and angle */
         *cphi = 1.0;                     /* to zero.  */
         *sphi = 0.0;
         return;
      } 
      else {                             /* Orientation change is small */
         *twitching = True;
         ca[1] = column_dot( T1, 3, T2, 2 );
         ca[2] = column_dot( T1, 1, T2, 3 );
         ca[3] = column_dot( T1, 2, T2, 1 );
         ca[4] = 0.0;
         d1 = dvmag( ca, 3 );
         if ( d1 > 0.00000000000001 )d1 = dvmnorm( ca, 3 );
         printf("Twitch\n");
         dvcopy( ca, 4, axis );
         *sphi = d1;
         *cphi = sqrt(1.0 - d1*d1);
         return;
      }
   }

   /*   Find angle   */

   dgetcolumn(T1, X, ca, 4);
   dgetcolumn(T1, Y, cb, 4);
   dgetcolumn(T2, X, da, 4);
   dgetcolumn(T2, Y, db, 4);
   if( dvdot(axis, 3, cb) > 0.9) { /*  rotation axis is close to j - use  */
      dcross4( axis, ca, cc );      /* change in i axis to find angle  */
      dcross4( axis, da, cd );    
   } 
   else {                          /* use  */
      dcross4( axis, cb, cc );      /* change in j axis to find angle  */
      dcross4( axis, db, cd );    
   }
   d1 = dvmag( cc, 3 );
   d2 = dvmag( cd, 3 );
   /*  fprintf( stdout, "d1, d2:  %11.6lf %11.6lf \n", d1, d2 ); */
   if ( d1 == 0.0 && d2 == 0.0 ) {
      *cphi = 1.0;
      *sphi = 0.0;
   } 
   else {
      dvsmy( cc, 3, 1.0/d1, cc );
      dvsmy( cd, 3, 1.0/d2, cd );   /* normalize */
      *cphi = dvdot(cc, 3, cd);     /* cos(phi) given by dot of two vectors */
      dcross4( cc, cd, dd );
      *sphi = dvmag(dd, 3);
      if( dvdot(dd, 3, axis) < 0.0) {
         *sphi = -(*sphi);
      }
   }
   return;
}

/* ---------------------------------------------------------------------------- */
/* Functions to rotate vector or frame about an arbitrary axis */

void RotateVector(double *v, double *axis, double sphi, double cphi, double *y)
{
   /* v= sin(angle)(axis X v) + Cos(angle)v + (1-Cos(angle))(v . axis) * axis
         See Robots for shearing sheep pg 128 
      Note that v and y can be the same 
   */

   double db[6],dd[6];
   double scal;

   validate_dvector( db, 1, 4);
   validate_dvector( dd, 1, 4);
   scal = (1.0 - cphi) * dvdot(axis, 3, v);
   dvsmy( axis, 3, scal, dd );  
   dcross4( axis, v, db);
   dvpiv( db, 3, sphi, dd, dd );
   dvpiv( v, 3, cphi, dd, y ); 
   y[4] = 0.0;
}

void RotateFrame(double **T1, double *axis, double sphi, double cphi,
     double **T2)
{
   /*  Routine to rotate coord frame T1 about axis by angle phi to give T2
       using quaternion rotation.  Note axis, T1 and T2 are all defined in
       terms of one (world) coord frame - axis is NOT defined in terms of T1.
       It is assumed that T1, and hence T2, are orthogonalized.  T1 can be
       the same array as T2.
       sphi, cphi are sin and cosine of 'phi'.
       Peter Kovesi    May 1985, adapted to C James Trevelyan July 1994
   */
   /*   Temporary local data    */
   double ca[6],cb[6],cc[6],cd[6]; /*   workspace vectors   */
   validate_dvector( ca, 1, 4);
   validate_dvector( cb, 1, 4);
   validate_dvector( cc, 1, 4);
   validate_dvector( cd, 1, 4);

#ifdef V_CHECK
   if ( !valid_dframe( T1 ) )
      nrerror("Invalid 1st frame: RotateFrame\n");
   if ( !valid_dframe( T2 ) )
      nrerror("Invalid 2nd frame: RotateFrame\n");
   if ( !valid_d4vector_b( axis ) )
      nrerror("Invalid axis vector: RotateFrame\n");
#endif

   dgetcolumn( T1, X, ca, 4 );
   RotateVector(ca, axis, sphi, cphi, cb);
   dgetcolumn( T1, Y, ca, 4 );
   RotateVector(ca, axis, sphi, cphi, cc);
   dcross4( cb, cc, cd); /* T2's k axis */
   dputcolumn( cb, 4, T2, X );
   dputcolumn( cc, 4, T2, Y );
   dputcolumn( cd, 4, T2, Z );
   dgetcolumn( T1, U, ca, 4 );
   dputcolumn( ca, 4, T2, U );
}

