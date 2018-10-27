/**
 This C program calculates the zeros of a Laguerre Polynomial
 L_n and their weight factors of the corresponding Gauss-Laguerre
 integration.
 Optionally it calculates the zeros of a Hermite Polynomial
 H_n and their weight factors of the corresponding Gauss-Hermite
 integration.

 The Laguerre Polynomial is
 L_n(x) = sum(i=0..n) (n choose i) (-x)^i/i!
 The zeros are defined to fullfil L_n(x_i)=0 for i=1,...,n.
 The weights are
 w_i = x_i/(n+1)^2/[L_(n+1)(x_i)]^2 = 1/x_i/[L_n'(x_i)]^2
 as given in section 25.4.45 in the "Handbook of Mathematical Functions"
 (after correcting an error of earlier printings). See also
 http://math.fullerton.edu/mathews/n2003/gausslaguerre/Gauss-LaguerreBib/Links/Gauss-LaguerreBib_lnk_2.html

 Usage from the UNIX shell after compilation:
   LagZero [-L | -H] n
 where n is a non-negative integer number denoting the order of the
 Laguerre or Hermite Polynomial.

 The output is a list of double column numbers similar to
 - the table 25.9 in the "Handbook  of Mathematical Functions", edited by
 M. Abramowitz and I. Stegun, which provides values up to n=15,
 - similar to the table in Appendix C of the book "Approximate calculate
 of integrals" by Vladimir Ivanoivich Krylov (McMillan, NY, 1962)
 which provides values up to n=32,
 - similar to the table on page 115 of my 1988 diploma thesis which covers
 the case n=50 (available from my publications web page).

 Example:
 unix>  LagZero -L 8
 1.7027963230510101e-01 3.6918858934163762e-01
 9.0370177679937991e-01 4.1878678081434306e-01
 2.2510866298661307e+00 1.7579498663717175e-01
 4.2667001702876588e+00 3.3343492261215170e-02
 7.0459054023934655e+00 2.7945362352257358e-03
 1.0758516010180996e+01 9.0765087733566494e-05
 1.5740678641278004e+01 8.4857467162725959e-07
 2.2863131736889265e+01 1.0480011748715021e-09

Compilation (Tested with gcc 4.1.1):

 gcc -std=gnu99 -O -o LagZero LagZero.c -lm

 The makefile entry would be:

LagZero: LagZero.c
        $(CC) -std=gnu99 -O -o $@ $< -lm

 The results might be compared against the output of running a Maple code
 of the following form:

 Digits := 16 ;
 n := 8 ; # for example
 fsolve(simplify(LaguerreL(n,0,x)),x,maxsols=n) ;
 
 Richard J. Mathar, http://www.strw.leidenuniv.nl/~mathar
**********************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* #define DEBUG */

/** Compute the ratio of the Laguerre Polynomial and its derivative
* The ratio of the value of the Laguerre Polynomial L_n^alpha and
* its derivative w.r.t. x, at the abscissa provided by the argument,
* is computed with the continued fraction representation of Shao et al
* [Math Comp 18 (88) (1964) 598, Equation (5.10)].
* @param[in] alpha the parameter of the generalized Laguerre Polynomial
* @param[in] n the index of the generalized Lagurerre Polynomial
* @param[in] x the abscissa point
* @return the ratio L_n^alpha(x)/L_n'^alpha(x)
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double LagRatio(const double alpha, const double x, const int n)
{
	int i=n ;
	double resul =0. ;
	while( --i >= 0)
		resul = (n-i+alpha)*(n-i)/(2*n-2*i-1+alpha-x-resul) ;
#ifdef DEBUG
	printf("ratio %lf at %lf %d\n",x/(n-resul),x,n ) ;
#endif
	return x/(n-resul) ;
}

/** Compute the ratio of the derivative of the Hermite Polynomial and its value.
* The inverse ratio of the value of the Hermite Polynomial H_n and
* its derivative w.r.t. x, at the abscissa provided by the argument,
* is computed with the continued fraction representation of
* eq (C3) in
* \latexonly
* \cite{MatharArxiv0306}.
* \endlatexonly
* \htmlonly
* <a href="http://arxiv.org/abs/math.NA/0306184">math.NA/0306184</a>.
* \endhtmlonly
* @param[in] n the index of the Hermite Polynomial
* @param[in] x the abscissa point
* @return the ratio H_n'(x)/H_n(x)
* @author Richard J. Mathar
* @since 2006-08-30
*/
static double HermRatio(const double x, const int n)
{
	int i=n ;
	double resul =0. ;
	while( --i >= 0)
		if ( i % 2 )
			resul = (n-i)/(2.*x-resul) ;
		else
			resul = (n-i)/(x-resul) ;
#ifdef DEBUG
	printf("ratio %lf at %lf %d\n",resul,x,n ) ;
#endif
	return resul ;
}

/** Compute a zero of the Laguerre Polynomial.
* A root of the Laguerre Polynomial L_n^alpha
* is computed with the third-order Newton method provided by Shao et al
* [Math Comp 18 (88) (1964) 598, Equation (5.2)], starting from a nearby
* estimate.
* @param[in] alpha the parameter of the generalized Laguerre Polynomial
* @param[in] n the index of the generalized Lagurerre Polynomial
* @param[in] x the initial guess of the root
* @return un update of the root fulfilling L_n^(alpha)(x)=0.
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double LagZeroNewton(const double alpha, double x, const int n)
{
	/* Loop until convergence to some double precision accuracy is reached */
	double xold ;
	for(int loop=0; loop<20 ;loop++)
	{
		/* Compute f/f' in the usual sense of the quadratic Newton method */
		const double ffprime= LagRatio(alpha,x,n) ;
		/* Use the value known up to here as a reference that determines
		* whether we have come close to a converged result.
		*/
		xold=x ;
		/* Use the cubically convergent Newton method for the update of x */
		x = xold -ffprime*(1.+0.5*(1.-(1.+alpha)/xold)*ffprime) ;
#ifdef DEBUG
		printf("%lf %lf %le\n",xold,x,x-xold) ;
#endif
		/* If double precision has been reached, leave the infinite loop and return */
		if ( x == xold)
			break ;
	}
	return x ;
}

/** Compute a zero of the Hermite Polynomial.
* A root of the Hermite Polynomial H_n
* is computed with the third-order Newton method provided by eq. (C1) of
* \latexonly
* \cite{MatharArxiv0306}.
* \endlatexonly
* \htmlonly
* <a href="http://arxiv.org/abs/math.NA/0306184">math.NA/0306184</a>.
* \endhtmlonly
* @param[in] n the index of the Hermite Polynomial
* @param[in] x the initial guess of the root
* @return un update of the root fulfilling H_n(x)=0.
* @author Richard J. Mathar
* @since 2006-08-30
*/
static double HermZeroNewton(double x, const int n)
{
	/* Loop until convergence to some double precision accuracy is reached */
	double xold ;
	for(int loop=0; loop<20 ;loop++)
	{
		/* Compute f'/f in the usual sense of the quadratic Newton method */
		const double fprimef= HermRatio(x,n) ;
#ifdef DEBUG
		printf("x=%lf n=%d loo=%d f'/f=%lf\n",x,n,loop,fprimef) ;
#endif
		/* Use the value known up to here as a reference that determines
		* whether we have come close to a converged result.
		*/
		xold=x ;
		/* Use the cubically convergent Newton method for the update of x */
		x = xold -(1.+xold/fprimef)/fprimef ;
#ifdef DEBUG
		printf("%lf %lf %le\n",xold,x,x-xold) ;
#endif
		/* If double precision has been reached, leave the infinite loop and return */
		if ( x == xold)
			break ;
	}
	return x ;
}

/** Return a guess of the i'th zero of the Laguerre Polynomial of order n.
* This uses Eq (5.6) in Shao et al, Math Comp 18 (1964) p 598.
* @param[in] n the index (order) of the polynomial
* @param[in] i the index of the zero, 0<=i<n .
* @return an estimate of the i'th root.
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double LagZeroGuess(const int n, const int i)
{
	/* The reference counts zeros as 1<=i<=n. So we shift
	* the definition by 1 for our use here: pi^2/n*[(4i+3)/8]^2 at alpha=0.
	*/
	return pow(M_PI*(0.5*i+0.375),2.)/n ;
}

/** Return a guess of the i'th non-negative zero of the Hermite Polynomial of order n.
* @param[in] n the index (order) of the polynomial
* @param[in] i the index of the zero, 0<=i<n .
* @return an estimate of the i'th root. This is the estimate of
*   eqs 22.15.3 and 22.15.4 of 
* \latexonly
* \cite{AS}.
* \endlatexonly
* \htmlonly
* the Handbook of Mathematical Function (edted by M Abramowitz, I Stegun).
* \endhtmlonly
* @author Richard J. Mathar
* @since 2006-08-30
*/
static double HermZeroGuess(const int n, const int i)
{
	if ( n %2)
		/* eq 22.15.4 of Abramowitz/Stegun solved for x */
		return 0.5*i*M_PI/sqrt(0.5*n) ;
	else
		/* eq 22.15.3 of Abramowitz/Stegun solved for x */
		return 0.5*(i+0.5)*M_PI/sqrt(0.5*n-0.5);
}

/** Lagrange extrapolation of f[0,1,2,3,4] to the next value.
* This is the five point interpolation formula, eq 25.2.15 in
* the "Handbook of Mathematical Functions" (eds: M. Abramowitz, I. Stegun).
* @param[in] f the five values in this vector are previous values at equidistant abscissa
* @return an extrapolation (estimate) to the next value in the vector,
*    equivalent to f[5].
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double LagrangeExtrap(const double *f)
{
	const double p=3. ;
	return (p+1.)*(p-1.)*p*(p-2.)/24.*f[0]
		-(p-1.)*p*(p-2.)*(p+2.)/6.*f[1]
		+(p+1.)*(p-1.)*(p-2.)*(p+2.)/4.*f[2]
		-(p+1.)*p*(p-2.)*(p+2.)/6.*f[3]
		+(p+1.)*(p-1.)*p*(p+2.)/24.*f[4] ;
}

/** Fill a vector with zeros of the Laguerre Polynomial.
* @param[in] n the index (order) of the polynomial
* @param[in,out] resul the storage for n results, each a zero
* @author Richard J. Mathar
* @since 2006-07-21
*/
static void LagZeroMain(const int n, double *resul)
{
	/* Eq (5.6) of the Shao reference provides initial guesses for the smallest zeros */
	for(int i=0; i < 2 ; i++)
		resul[i] = LagZeroGuess(n,i) ;
#ifdef DEBUG
	/* printf("%lf %lf estimated zeros from a1=%lf a0 = %lf\n",resul[0],resul[1],a1,a0) ; */
	printf("%lf %lf estimated zeros\n",resul[0],resul[1]) ;
#endif

	/* Improve both estimates with the Newon Method */
	for(int i=0 ; i < 2 ; i++)
		resul[i] = LagZeroNewton(0.,resul[i],n) ;
#ifdef DEBUG
	printf("%lf %lf zeros\n",resul[0],resul[1]) ;
#endif

	resul[2] = LagZeroGuess(n,2) ;
	resul[2] = LagZeroNewton(0.,resul[2],n) ;
#ifdef DEBUG
	printf("%lf zero\n",resul[2]) ;
#endif

	/* The rest of the zeros are hopefully almost equidistant such that we can
	* bootstrap them by extrapolation from those already found.
	*/
	for(int i=3; i < n; i++)
	{
		if ( i <= 4 )
			resul[i] = LagZeroGuess(n,i) ;
		else
			resul[i] = LagrangeExtrap(&resul[i-5]) ;
		resul[i] = LagZeroNewton(0.,resul[i],n) ;
#ifdef DEBUG
		printf("%lf zero\n",resul[i]) ;
#endif
	}
}

/** Fill a vector with zeros of the Hermite Polynomial.
* @param[in] n the index (order) of the polynomial
* @param[in,out] resul the storage for (n+1)/2 results, each a zero
* @author Richard J. Mathar
* @since 2006-08-30
*/
static void HermZeroMain(const int n, double *resul)
{
	int i=0 ;
	for( ; i < (n+1)/2 && i <= 2; i++)
	{
		resul[i] = HermZeroGuess(n,i) ;
		/* Improve both estimates with the Newton Method */
		resul[i] = HermZeroNewton(resul[i],n) ;
	}

	/* The rest of the zeros are hopefully almost equidistant such that we can
	* bootstrap them by extrapolation from those already found.
	*/
	for( ; i < (n+1)/2; i++)
	{
		if ( i <= 4 )
			resul[i] = HermZeroGuess(n,i) ;
		else
			resul[i] = LagrangeExtrap(&resul[i-5]) ;
		resul[i] = HermZeroNewton(resul[i],n) ;
#ifdef DEBUG
		printf("%lf zero\n",resul[i]) ;
#endif
	}
}

/** Compute a vector with zeros of the Laguerre Polynomial.
* @param n the index (order) of the polynomial
* @return a vector of length n with the zeros sorted in increasing order.
*   This is actually a pointer to the first value and should be free()'d by the
*   caller after having used the result.
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double *LagZero(const int n)
{
	double * resul = NULL ;
	if ( n < 0 )
	{
		fprintf(stderr,"Invalid argument n = %d in "__FILE__" line %d\n",n,__LINE__) ;
		return NULL ;
	}
	resul = (double*)malloc(n*sizeof(double)) ;
	switch(n)
	{
	case 0:
		/* L_0(x)=1 which has no zeros. */
		break ;
	case 1 :
		/* L_1(x)=1-x which has the zero x=1, obviously */
		resul[0]=1.0 ;
		break ;
	case 2:
		/* L_2(x)=1-2x+x^2/2, which has the zeros x=2+-sqrt(2) */
		resul[0]=2.-M_SQRT2 ;
		resul[1]=2.+M_SQRT2 ;
		break ;
	default:
		LagZeroMain(n,resul) ;
	}
	return resul ;
}

/** Compute a vector with zeros of the Hermite Polynomial.
* @param n the index (order) of the polynomial
* @return a vector of length (n+1)/2 with the non-negativer zeros sorted in increasing order.
*   This is actually a pointer to the first value and should be free()'d by the
*   caller after having used the result.
* @author Richard J. Mathar
* @since 2006-08-30
*/
static double *HermZero(const int n)
{
	double * resul = NULL ;
	if ( n < 0 )
	{
		fprintf(stderr,"Invalid argument n = %d in "__FILE__" line %d\n",n,__LINE__) ;
		return NULL ;
	}
	resul = (double*)malloc((1+n)/2*sizeof(double)) ;
	switch(n)
	{
	case 0:
		/* H_0(x)=1 which has no zeros. */
		break ;
	case 1 :
		/* H_1(x)=2x which has the zero x=0, obviously */
		resul[0]=0.0 ;
		break ;
	case 2:
		/* H_2(x)=4x^2-2, which has the zeros x=+-1/sqrt(2) */
		resul[0]= M_SQRT1_2 ;
		break ;
	default:
		HermZeroMain(n,resul) ;
	}
	return resul ;
}

/** Compute derivative of Laguerre Polynomial
* @param n the order of the polynomial
* @param[in] x the abscissa at which the derivative is requested
* @return the value of L_n'(x)
* @author Richard J. Mathar
* @since 2006-07-21
*/
static double LagDeriv(const int n, const double x)
{
	/* The derivative is L_n'(x)=sum(j=0..n-1) (n choose j+1) (-)^(j+1)/j! z^j
	*/
	double resul = -n,
		fac = resul ;
	switch(n)
	{
		/* for n=0 and 1, the derivatives are 0 and -1 and already correct
		* from the setting of 'resul' above
		*/
	case 0:
	case 1:
		break ;
	default:
		for(int j=1; j < n; j++)
		{
			fac *= -x*(double)(n-j)/(j*(j+1)) ;
			resul += fac ;
		}
	}
	return resul ;
}

/** Compute the value of the Hermite Polynomial
* @param n the order of the polynomial
* @param[in] x the abscissa at which the value is requested
* @return the value of H_n(x)
* @author Richard J. Mathar
* @since 2006-08-30
*/
static double Herm(const int n, const double x)
{
	/* The result is H_n(x)= n! sum(j=0..n/2) (-)^j/j!/(n-2j)! (2x)^(n-2j)
	*/
	double resul = pow(2.*x,(double)n) ,
		fac = resul ;
	if ( x == 0.)
	{
		/* case 22.4.8 of Abramowitz/Stegun */
		const int m=n/2 ;
		switch (n % 2)
		{
		case 0 :
			resul= 1. ;
			for(int j= 1+m; j <= n ; j++)
				resul *= j;
			return (m %2 ) ? -resul : resul ;
		case 1 :
			return 0. ;
		}
	}
	switch(n)
	{
		/* for n=0 and 1, the values are 1 and 2x.
		* from the setting of 'resul' above
		*/
	case 0:
		resul = 1. ;
		break ;
	default:
		for(int j=1; 2*j <= n; j++)
		{
			fac *= -0.25*(n+2-2*j)*(n+1-2*j)/(x*x*j) ;
			resul += fac ;
		}
	}
#ifdef DEBUG
	printf("H(n=%d x=%lf)= %lf\n",n,x,resul) ;
#endif
	return resul ;
}

/** The main program from the UNIX shell.
* @param[in] argc
* @param[in] argv the two command line arguments. The first is a mandatory
*   switch indicative of whether we are handling the Laguerre (-L) or Hermite (-H)
*   polynmials. The second is the polynomial order.
* @return 1 if command line parameters or call were in error, 0 otherwise.
* @author Richard J. Mathar
* @since 2006-07-21
*/
int main(int argc, char* argv[])
{
	/** The variable \c n is the order of the polynomial. 
	*/
	int n= -1 ;
	/** The variable \c polytpe is the kind of the polynomial. 0 for the Laguerre case
	* 1 for the Hermite case.
	*/
	int poltype= -1 ;
	double *x = NULL,
	  *w = NULL ;

	/* The command line parameter
	*/
	char oc ;

	/* Check usage: need one command line parameter, which is
	* the order of the polynomial, n>=0 .
	*/
	if ( argc < 3 )
	{
		fprintf(stderr,"usage: %s [-L | -H] n\n",argv[0]) ;
		return 1 ;
	}

        while (  (oc=getopt(argc,argv,"LH")) != -1 )
        {
                switch(oc)
                {
                case 'L' :
			poltype=0 ;
                        break ;
                case 'H' :
			poltype=1 ;
                        break ;
                case '?' :
                        fprintf(stderr,"Invalid command line option %c\n",optopt) ;
                        break ;
                }
        }

	n = atoi(argv[optind]) ;
	if ( n < 0 )
	{
		fprintf(stderr,"Usage error need non-negative n, got %d\n",n) ;
		return 1 ;
	}

#ifdef DEBUG
	printf("n= %d\n",n) ;
#endif
	switch( poltype)
	{
	case 0:
		/* Obtain the vector of the zeros by calling LagZero().
		*/
		x = LagZero(n) ;

		/* provide storage for the weights of the Gauss-Laguerre integration
		*/
		w = (double*)malloc(n*sizeof(double)) ;
		/* Compute the weights
		*/
		for(int i=0; i < n ; i++)
		{
			const double Lprime = LagDeriv(n,x[i]) ;
#ifdef DEBUG
			printf("%d %lf\n",i,Lprime) ;
#endif
			w[i] = 1./(Lprime*Lprime*x[i]) ;
		}

		/* Print result to stdout
		*/
		for(int i=0 ; i < n ; i++)
			printf("%.16le %.16le\n",x[i],w[i]) ;

		break ;
	case 1:
		/* Obtain the vector of the zeros by calling HermZero().
		*/
		x = HermZero(n) ;

		/* provide storage for the weights of the Gauss-Laguerre integration
		*/
		w = (double*)malloc((1+n)/2*sizeof(double)) ;
		/** Compute the weights with eq (C4) of
		* \latexonly
		* \cite{MatharArxiv0306}.
		* \endlatexonly
		* \htmlonly
		* <a href="http://arxiv.org/abs/math.NA/0306184">math.NA/0306184</a>.
		* \endhtmlonly
		*/
		for(int i=0; i < (1+n)/2 ; i++)
		{
			/* H contains n*H(n-1,x) */
			const double H = n*Herm(n-1,x[i]) ;
#ifdef DEBUG
			printf("i=%d n-1=%d x=%lf %lf\n",i,n-1,x[i],H) ;
#endif
			/* start with 2^-1*sqrt(pi)/n^2/H(n-1,x)^2 */
			w[i] = 1./(M_2_SQRTPI*H*H) ;
			for (int j=1; j <=n ; j++)
				w[i] *= 2*j ;
		}
		/* Print result to stdout
		*/
		for(int i=0 ; i < (n+1)/2 ; i++)
			printf("%.16le %.16le\n",x[i],w[i]) ;
		break ;
	default:
		fprintf(stderr,"Invalid choice of polynomial type\n") ;
		return 1 ;
	}

	/* Some cleanup (which is not stricly needed but shown to demonstrate
	* the interfaces).
	*/
	free(x) ;
	free(w) ;
	return 0 ;
}
