/* This is a sample problem taken from Moshier's Cephes library which
 * includes a port of LMDIF. The quantities MACHEP and DWARF would
 * normally be replaced with the standard IEEEE machine constants
 * for macheps and min_float.
 */
#include <math.h>
#include <stdio.h>
#include "cminpak.h"
#include "dpmpar.h"

#define N 4
#define M 4
double tol = 1.0e-14;
double epsfcn = 1.0e-15;
double factor = 0.1;
double x[N] = {0.0};
double fvec[M] = {0.0};
int msk[N] = {1,1,1,1};
int ipvt[N] = {0};
int maxfev = 200 * (N+1);
/* resolution of arithmetic */
double MACHEP = 1.2e-16;
/* smallest nonzero number */
double DWARF = 1.0e-38;
void fcn(int m,int n,double x[],double fvec[],int *iflag);

int main()
{
int m,n,info;
int mode,nfev;
int iflag, ecode;

n = N;
m = M;
printf( "initial x\n" );
printf("%lf  %lf  %lf  %lf\n",x[0],x[1],x[2],x[3]);
fcn(m,n,x,fvec,&iflag);
printf( "initial function\n" );
printf("%lf  %lf  %lf  %lf\n",fvec[0],fvec[1],fvec[2],fvec[3]);
/* Call lmdif. */
mode = 1;
info = 0;
nfev = 0;

ecode=lmdif0(fcn,m,n,x,msk,fvec,tol,&info,&nfev);

printf( "%d function evaluations\n", nfev );
/* display solution and function vector */
printf( "x\n" );
printf("%lf  %lf  %lf  %lf\n",x[0],x[1],x[2],x[3]);
printf( "fvec\n" );
printf("%lg  %lg  %lg  %lg\n",fvec[0],fvec[1],fvec[2],fvec[3]);
printf( "function norm = %.15e\n", enorm(m, fvec) );
return 0;
}

/***********Sample of user supplied function****************
 * m = number of functions
 * n = number of variables
 * x = vector of function arguments
 * fvec = vector of function values
 * iflag = error return variable
 */
void fcn(int m,int n,double x[],double fvec[],int *iflag)
{
double temp;

/* an arbitrary test function: */
fvec[0] = 1.0 - 0.3 * x[0] + 0.9 * x[1] - 1.7 * x[2] +log(1.5+x[3]);
fvec[1] = sin(-4.0*x[0] ) - 3.0 * x[1] + 0.1 * x[2] + x[3]*x[3];
temp = (x[2] + 2.0) * x[2] * x[1];
fvec[2] = 0.5 * x[1] - sin( x[2] + 1.0 ) + temp + .3 * x[3];
fvec[3] = x[0]*x[1] + x[1]*x[2] + x[0]*x[2] - x[3]*x[3];
}
