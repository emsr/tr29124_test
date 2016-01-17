#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

			/* define HAVE_NO_CBRT if your math library does not contain a double cbrt(double)
			function which returns the cube root of a single real number, returning a negative
			number for negative arguments */
/*#define HAVE_NO_CBRT */
			/* define HAVE_NO_SINCOS if your math library does not contain a void sincos(double,double *,double*)
			function which returns the sine and cosine a single real number at the same time */
/*#define HAVE_NO_SINCOS */
/**********************************************************************
 Solve a cubic equation with real-valued coefficients a[0..2]
 a[0]+a[1]*z+a[2]*z^2+z^3=0
 for 'z' (returned in z[0..2]) and return the number of real-valued results
 (1 or 3).  z[][0] contain the real parts of the results, z[][1] the imaginary parts.

 Method: see formulas in Abramowitz/Stegun chapter 3.82
         http://mathworld.wolfram.com/CubicEquation.html  eqs (43)-(45)

 Richard J Mathar, 11 Oct 2000
*************************************************************************/
#define M_SQRT3 1.73205080756887729352744634151
#define M_SQRT3_2 .866025403784438646763723170755
int cuberoots(const double a[3], double z[3][2])
{
	const double a2third=a[2]/3. ,
	             q=a[1]/3.-a2third*a2third ;
	double r=(a[1]*a2third-a[0])/2.-pow(a2third,3.),	/* check that pow(negative_here,3.) delivers correct sign .... */
	       disc=pow(q,3.)+r*r ;			/* discriminant */
					/* consider the degenerate case of a[0]=0 where it all becomes a quadratic eqn */
	if( a[0] == 0.)
	{
		const double a2half=a[2]/2. ;
		z[0][0]=z[0][1]=0. ;
		disc= a2half*a2half-a[1] ;
		if( disc >= 0.)			/* two remaining real roots of quadratic equation */
		{
			z[1][0]= -a2half+sqrt(disc) ;
			z[2][0] = -a[2]-z[1][0] ;
			z[1][1]=z[2][1]=0. ;
			return 3 ;
		}
		else				/* two remaining complex conjugate roots of quadratic equation */
		{
			z[1][0]= z[2][0]= -a2half ;
			z[2][1]= - ( z[1][1] = sqrt(-disc) ) ;
			return 1 ;
		}
	}

	if (disc >0.)				/* one real and a pair of complex conjugate roots */
	{
		double s[2] ;
		disc=r+sqrt(disc) ;
						/* better not use s[1]=cbrt(r-sqrt(disc)) because the branches of the
						cube roots must be chosen such that s[0]*s[1]=-q */
		s[1]= -q/(s[0]=cbrt(disc)) ;
		z[0][0]=s[0]+s[1]-a2third ;
		z[0][1]=0. ;
		z[1][0]=z[2][0]= -(s[0]+s[1])/2.-a2third ;
		z[2][1]= -(z[1][1]= M_SQRT3_2*(s[0]-s[1]) ) ;
		return 1 ;
	}
	else				/* three real valued roots */
	{
		struct { double re; double im; double modul ; double phase ; } s[2] ; /* two complex numbers */
		s[0].re = r ;
		s[0].im = sqrt(-disc) ;		/* s[0]=r+(disc)^(1/2) with disc <0 */
		s[0].phase= atan2(s[0].im,s[0].re)/3. ;	/* compute s[0]^(1/3) in spherical/Euler representation */
		s[0].modul= pow(r*r-disc,1./6.) ;	/* modulus of s[0] is sqrt(re*re+im*im). compute modulus
							of s[0]^(1/3) =(s[0].re^2+s[0].im^2)^(1/6)  where s[0].im^2= -disc */
							/* and back to Cartesian coordinates */
		sincos(s[0].phase,&s[0].im,&s[0].re) ;
		s[0].re *= s[0].modul ;
		s[0].im *= s[0].modul ;
							/* compute s[1]= -q/s[0] */
		s[1].re = -q*s[0].re/(s[0].modul*s[0].modul) ;
		s[1].im = q*s[0].im/(s[0].modul*s[0].modul) ;

							/* z[0]= s[0]+s[1] -a[2]/3 */
		z[0][0]= s[0].re+s[1].re-a2third ;
							/* z[1]= -{s[0]+s[1]+i sqrt(3) (s[0]-s[1])}/2 -a[2]/3 with
							i*(s[0]-s[1])= -(Im s[0]-Im s[1]); only real parts needed.... */
		z[1][0]= -(s[0].re+s[1].re+M_SQRT3*(s[0].im-s[1].im))/2.-a2third ;
							/* z[2]= -{s[0]+s[1]-i sqrt(3) (s[0]-s[1])}/2 -a[2]/3 */
		z[2][0] = -a[2]-z[0][0]-z[1][0] ;
		z[0][1]=z[1][1]=z[2][1] =0. ;
		return 3 ;
	}
}
#undef M_SQRT3_2
#undef M_SQRT3

#ifdef HAVE_NO_CBRT
/******************************************************************
 Implement a local cbrt() function (cube root) if not provided by the local library.
 Returns "sign-extended" negative numbers if the argument is negative.
 Not needed for Sun's or HP's native compilers, for example.
 Richard J. Mathar, 11 Oct 2000
******************************************************************/
static double cbrt(const double x)
{
	if( x < 0.)
		return -pow(-x,1./3.) ;
	else
		return pow(x,1./3.) ;
}
#endif /* HAVE_NO_CBRT */

#ifdef HAVE_NO_SINCOS
/***************************************************************************
 Auxiliary replacement for the sincos() function that might be missing
 in some mathematical environments.
 Available with   -lsunmath -lm with the Sun compilers.
 Richard J. Mathar, 16 Oct 2000
***************************************************************************/
static void sincos(const double x, double *s, double *c)
{
	*s=sin(x) ;
	*c=cos(x) ;
}
#endif /* HAVE_NO_SINCOS */
