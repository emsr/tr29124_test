/*******************************************************************
compute the coefficients a[i], i=0,..n, of the Jacobi polynomial
P(1-2x) [order n, upper parameters (alpha,beta)=(1,0)]
P(1-2x)=(n+1)*sumover i from 0 of( (-n)sub i *(n+2)sub i*x^i/[i!(i+1)!]
with (a)sub i denoting Pochhammer's symbol, '^' the exponentiation and
'!' the factorial.

Literature: M. Abramowitz, I. Stegun (edts), 'Handbook of Mathematical
Functions with Formulas, Graphs and Mathematical TAbles', in the series
'Dover books in Advanced mathematics', Dover Publications, INC, New York,
9th Dover printing, 1973, ISBN 0-486-612672-4,
formulas (22.5.42),(15.4.1,(6.1.22)
The result must be the 1st moment (k=1) entrances of Table 25.8 therein.

usage:
	gauss n

As a function of the NAG-library is used, compile it on a Convex machine with
	cc -c gauss.c
	fc -o gauss gauss.o -lnag

R.J. Mathar, FZR/FIA, 15.06.92
********************************************************************/
		/* maximum order */
#define NMAX 50
		/* to declare 'stderr' */
#include <stdio.h>
		/* to declare 'exit() */
#include <stdlib.h>
main(int argc, char *argv[])
{
	int i,j,
	    scale=1,			/* FOTRAN .true. hopefully */
	    n,
	    ifail=0;
	double zero[NMAX][2],	/* zeroes of the Jacobi polynomial */
	       weight[NMAX],	
	       allcffs[NMAX+1][NMAX+1],
	       worksp[2*NMAX+2];
	extern void c02agf_(double *,int *,int *,double (*)[2],double *,int *) ;
	void Coffs(const, double*) ;

	sscanf(argv[1],"%d",&n) ;
	if( n > NMAX)
	{
		fprintf(stderr,"n greater than %d: recompile\n",n) ;
		exit(1) ;
	}
	for(j=0;j<=n;j++)
		Coffs(j,allcffs[j]) ;
			/* note: we've computed the coefficients in the backward
			order defined and needed in the c02agf_() function. As we
			don't need allcffs[n][0..n] any more, we can
			reverse the order here on site. */
	for(i=0;i<=n/2;i++)
	{
		register double tmp=allcffs[n][i];
		allcffs[n][i]=allcffs[n][n-i] ;
		allcffs[n][n-i]=tmp;
	}
			/* compute all zeros if (sum i=0 to n) allcffs[n][i]x^i
			-> zero[0,...,n-1] */
	c02agf_(allcffs[n],&n,&scale,zero,worksp,&ifail) ;
	for(i=0;i<n;i++)
	{
		weight[i]=0. ;
		for(j=0;j<n;j++)
		{
			double horner(int,double,double*) ;
			register double P10j=horner(j,zero[i][0],allcffs[j]) ;
			weight[i] += P10j*P10j*(1+j) ;
		}
		weight[i]= .5/weight[i] ;
	}
			/* revers order of output .. */
	printf("\t[%d]={\n\t",n) ;
	for(i=n-1;i>0;i--)
		printf("%.15le,\n\t",zero[i][0]) ;
	printf("%15le};\n",zero[0][0]) ;
	printf("\t[%d]={\n\t",n) ;
	for(i=n-1;i>0;i--)
		printf("%.15le,\n\t",weight[i]) ;
	printf("%15le};\n",weight[0]) ;
}
/*************************************************
 compute the polynomial (sum i=0up to n)
 coeff[i]*x^i with help of the Horner scheme
 R.J. Mathar, FZR/FIA, 15.06.92
*****************************************************/
double horner(int n, double x, double Coffs[])
{
	int rec;
	double sum=Coffs[n];
	for(rec=n-1;rec>=0;rec--)
		sum=sum*x+Coffs[rec];
	return sum;
}
/**************************************************
 Compute the coefficients coeff[0,...,n] for the
 Jacobi polynomial P(1,0)n of order 'n'.
 P(alpha=1,beta=0,n)(1-2x)=(sum i=0 up to n) coeff[i]*x^i.
*************************************************/
void Coffs(const n, double coeff[])
{
	double nom,		/* Pochhammer's symbol (-n)sub i*(n+2)sub i, numerator */
	       den=1.;		/* factorial i!, denominator */
	int i;
	coeff[0]=nom=n+1 ;
	for(i=1;i<=n;i++)
	{
		nom *= (i-n-1)*(n+i+1) ;
		den *= i*(i+1) ;
		coeff[i]=nom/(double)den;
	}
}
	
