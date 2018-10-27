/* deint.c -- integration of function by double exponential transform
 *  from Applied Numerical Methods with Software, by Nakamura.
 *
 *  (C) 2001, C. Bond. All rights reserved.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>

double deint(double (*f)(double),int n,double a,double b)
{
    double h,ss,pi,z,exz,hcos,hsin,s,dxdz,x,w;
	int k;

	n = n/2;
    h = 5.0/n;    /* 5.0 is rough limit of K in exp(exp(K)) */

    pi = M_PI_4;
	ss = 0.0;
	for (k=-n;k<=n;k++) {
		z = h * (double)k;
		exz = exp(z);
		hcos = exz+1.0/exz;
		hsin = exz-1.0/exz;
        s = exp(pi*hsin);
		w = s + 1.0/s;
		x = (b*s+a/s)/w;     /* transformed abscissa  */
		if (x != a && x != b) {
            dxdz = 2*(b-a)*pi*hcos/(w*w); /* transformed weight */
			ss += h * f(x)*dxdz;
		}
	}
	return ss;
}

double deint2(double (*f)(double),int n,double a,double b)
{
    double h,ss,pi,z,exz,hcos,hsin,s,dxdz,x1,x2,w;
	int k;

	n = n/2;
    h = 5.0/n;    /* 5.0 is rough limit of K in exp(exp(K)) */

    pi = M_PI_4;
	ss = 0.5*f((a+b)/2);
	for (k=-n;k<0;k++) {
		z = h*(double)k;
		exz = exp(z);
		hcos = exz+1.0/exz;
		hsin = exz-1.0/exz;
        s = exp(pi*hsin);
		w = s + 1.0/s;
		dxdz = hcos/(w*w);      /* transformed weight   */  
		x1 = (b*s+a/s)/w;       /* transformed abscissa */
		x2 = (a*s+b/s)/w;       /*     "         "      */
		if (x1 != a && x1 != b) 
			ss += dxdz * f(x1);
		if (x2 != a && x2 != b)
			ss += dxdz * f(x2);
	 }
    return 2*(b-a)*pi*h*ss;
}

double deint3(double (*f)(double),int n,double a,double b)
{
    double h,ss,pi,z,exz,hcos,hsin,s,dxdz,x,w;
	int k;

	n = n/2;
    h = 5.0/n;    /* 5.0 is rough limit of K in exp(exp(K)) */

    pi = M_PI_4;
	ss = 0.0;
	for (k=-n;k<=n;k++) {
		z = h * (double)k;
		exz = exp(z);
		hcos = exz+1.0/exz;
		hsin = exz-1.0/exz;
        s = exp(pi*hsin);
		w = s + 1.0/s;
		x = (b*s+a/s)/w;        /* transformed abscissa  */
		if (x != a && x != b) {
			dxdz = hcos/(w*w);  /* transformed weight    */
			ss += f(x) * dxdz; 
		}
	}
    printf ("1st subtotal: %15.15lf\n", h*2*(b-a)*pi*ss);

	for (k=-n;k<=n;k++) {
		z = h * ((double)(k)+0.5);      /* use in-between points this pass */
		exz = exp(z);
		hcos = exz+1.0/exz;
		hsin = exz-1.0/exz;
        s = exp(pi*hsin);
		w = s + 1.0/s;
		x = (b*s+a/s)/w;        /* transformed abscissa  */
		if (x != a && x != b) {
			dxdz = hcos/(w*w);  /* transformed weight    */
			ss += f(x) * dxdz; 
		}
	}
     printf ("2nd subtotal: %15.15lf\n", h*(b-a)*pi*ss);

    return h*(b-a)*pi*ss;
}
/* Progressive version */
double deint4(double (*f)(double),int n,double a,double b)
{
    double h,ss,ss1,ss2,pi,z,exz,hcos,hsin,s,dxdz,x1,x2,w;
	int k,l;

	n = n/2;
    h = 5.0/n;    /* 5.0 is rough limit of K in exp(exp(K)) */

    pi = M_PI_4;
	ss = 0.5*f((a+b)/2);
	ss1 = 0.0;
	ss2 = 0.0;
	for (k=-n;k<0;k++) {
		z = h*(double)k;
		exz = exp(z);
		hcos = exz+1.0/exz;
		hsin = exz-1.0/exz;
        s = exp(pi*hsin);
		w = s + 1.0/s;
		dxdz = hcos/(w*w);      /* transformed weight   */  
		x1 = (b*s+a/s)/w;       /* transformed abscissa */
		x2 = (a*s+b/s)/w;       /*     "         "      */
		if (x1 != a && x1 != b) 
			ss1 += dxdz * f(x1);
		if (x2 != a && x2 != b)
			ss2 += dxdz * f(x2);
	}
    printf("1, %15.15lf, %15.15lf\n",2*(b-a)*pi*h*ss1,2*(b-a)*pi*h*ss2);

	for (l=0;l<3;l++) {

		for (k=-n;k<0;k++) {
			z = h*((double)k+0.5);
			exz = exp(z);
			hcos = exz+1.0/exz;
			hsin = exz-1.0/exz;
            s = exp(pi*hsin);
			w = s + 1.0/s;
			dxdz = hcos/(w*w);      /* transformed weight   */ 
			x1 = (b*s+a/s)/w;       /* transformed abscissa */
			x2 = (a*s+b/s)/w;       /*     "         "      */
			if (x1 != a && x1 != b) 
				ss1 += dxdz * f(x1);
			if (x2 != a && x2 != b)
				ss2 += dxdz * f(x2);
		}
        printf("%d, %15.15lf, %15.15lf\n",l+2,(b-a)*pi*h*ss1,(b-a)*pi*h*ss2);
		n *= 2;
		h /= 2;
	}
    return 2*(b-a)*pi*h*(ss+ss1+ss2);
}

