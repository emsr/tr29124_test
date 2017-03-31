/* bairstow.c -- Bairstow's complex root finder
 *
 *  (C) 1991, C. Bond
 *
 *  Finds all real and complex roots of polynomials
 *  with real coefficients.
 *
 *  Features:
 *      o Global, self-starting method,
 *      o Does not require initial estimates,
 *      o Finds all roots and quadratic factors,
 *
 *  Caveats:
 *      All polynomial root finders have problems
 *      maintaining accuracy in the presence of
 *      repeated roots.
 *      This is no exception!
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MAX_TERMS 21
int precision_error_flag;
double a[MAX_TERMS],b[MAX_TERMS],c[MAX_TERMS],d[MAX_TERMS];

void find_poly_roots(int n)
{
		double r,s,dn,dr,ds,drn,dsn,eps;
	int i,iter;

	r = s = 0;
	dr = 1.0;
	ds = 0;
	eps = 1e-14;
	iter = 1;

	while ((fabs(dr)+fabs(ds)) > eps) {
		if ((iter % 200) == 0) {
			r=(double)rand()/16000.;
		}
		if ((iter % 500) == 0) {
			eps*=10.0;
			precision_error_flag=1;
			printf("Loss of precision\n");
		}
		b[1] = a[1] - r;
		c[1] = b[1] - r;

		for (i=2;i<=n;i++){
			b[i] = a[i] - r * b[i-1] - s * b[i-2];
			c[i] = b[i] - r * c[i-1] - s * c[i-2];
		}
		dn=c[n-1] * c[n-3] - c[n-2] * c[n-2];
		drn=b[n] * c[n-3] - b[n-1] * c[n-2];
		dsn=b[n-1] * c[n-1] - b[n] * c[n-2];

		if (fabs(dn) < 1e-16) {
			dn = 1;
			drn = 1;
			dsn = 1;
		}
		dr = drn / dn;
		ds = dsn / dn;

		r += dr;
		s += ds;
		iter++;
	}
	for (i=0;i<n-1;i++) 
		a[i] = b[i];
	a[n] = s;
	a[n-1] = r;
}


void main()
{
	int i,n,order=0;
	double tmp;

	while ((order < 2) || (order > MAX_TERMS-1)) {
		printf("Polynomial order (2-20): ");
		scanf("%d",&order);
	}
	printf("Enter coefficients, high order to low order.\n");
	for (i=0;i<=order;i++){
		printf("a(%d)=",i);
		scanf("%lf",&a[i]);
		if (i==0) {
			tmp=a[i];
			a[i]=1.0;
		}
		else a[i]/=tmp;
		d[i]=a[i];
	}
	b[0]=c[0]=1.0;
	n=order;
	precision_error_flag=0;
	while (n > 2) {
		find_poly_roots(n);
				n-=2;
	}
	printf("The quadratic factors are:\n");
	
	for (i=order;i>=2;i-=2) /* print quadratics */
		printf("t^2  %+.15lg t  %+.15lg\n",a[i-1],a[i]);
	if ((n % 2) == 1)
		printf("The linear term is: \nt  %+.15lg\n",a[1]);
}

