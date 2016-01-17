#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdarg>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "poly3D.h"

/* availability of the factorials etc
*/
#include <gsl/gsl_sf.h>

/* availability of GSL_IS_ODD() for example.
*/
#include <gsl/gsl_math.h>

using namespace std ;

#ifdef INCLUDE_SUPERFLUOUS

/**Default creator with given 3D coordinates
*/
point3D::point3D( const double cart[TURB3D_DIM], double radius =1.) : rad(radius)
{
	memcpy(xyz,cart,sizeof(double[TURB3D_DIM]) ) ;
}
#endif

/** Default Ctor at the origin of coordinates
*/
point3D::point3D( )
{
	xyz[0]=xyz[1]=xyz[2]=0. ;
}

/** Ctor from individual cartesian coordinates
*/
point3D::point3D( const double X, const double Y, const double Z)
{
	xyz[0] = X ;
	xyz[1] = Y ;
	xyz[2] = Z ;
}

/** Ctor with given 3D coordinates
*/
point3D::point3D( const double cart[TURB3D_DIM])
{
	memcpy(xyz,cart,sizeof(double[TURB3D_DIM]) ) ;
}

/** Copy constructor.
*/
point3D::point3D( const point3D & orig)
{
	memcpy(xyz,orig.xyz,sizeof(double[TURB3D_DIM]) ) ;
}

/** Distance to the origin
*/
double point3D::dist() const
{
	/* distance to polar axis along the equatorial plane
	*/
	const double equat = hypot(xyz[0],xyz[1]) ;

	/* radial distance to origin 
	*/
	return hypot(equat,xyz[2]) ;
}

/** Distance to another point
*/
double point3D::dist( const point3D &oth) const
{
	/* distance to polar axis along the equatorial plane
	*/
	const double equat = hypot(xyz[0] - oth.xyz[0],xyz[1] -oth.xyz[1]) ;

	/* radial distance to origin 
	*/
	return hypot(equat,xyz[2] -oth.xyz[2]) ;
}

/** Scale the points such that those with distance 'radius'
* appear to be at distance 1
*/
void point3D::scale(const double radius)
{
	for(int i=0; i < TURB3D_DIM; i++)
		xyz[i] /= radius ;
}

/** Scale components so the new length equals a given number
* @param newleng the distance to the origin (or vector length) on return
*/
void point3D::normalize(const double newleng)
{
	const double len = dist() ;
	scale(len/newleng) ;
}

/** Convert point to spherical angular coordinates.
* @param rtp radius, theta and phi on return.
*       r = sqrt(x^2+y^2+z^2) ;
* 	x=r sin(theta) cos(phi)
* 	y=r sin(theta) sin(phi)
* 	z=r cos(theta)
*/
void point3D::phithet( double rtp[3]) const
{
	/* distance to polar axis along the equatorial plane
	*/
	const double equat = hypot(xyz[0],xyz[1]) ;

	/* radial distance to origin 
	*/
	rtp[0] = hypot(equat,xyz[2]) ;

	/* azimuthal angle phi in the range 0 to 2 pi
	*/
	rtp[2] = atan2(xyz[1],xyz[0]) ;

	/* polar angle theta in the range 0 to pi
	* equat=r sin(theta)
	*/
	rtp[1] = atan2(equat,xyz[2]) ;
}

/** Assigment Operator.
* @param right the right hand side of the assignment
* @return the new vector that results
*/
point3D & point3D::operator= (const point3D & right)
{
	if ( this != & right)
	{
		memcpy(xyz,right.xyz,sizeof(double[TURB3D_DIM]) ) ;
	}
	return *this ;
}

/** Multiply by a constant in the sense that this is a vector shortened/lengthed.
* @param c the constant to multply with
* @return the new vector that results
*/
point3D & point3D::operator*= (const double c)
{
	for(int d=0; d < TURB3D_DIM ; d++)
		xyz[d] *= c ;
	return *this ;
}

/** Add a vector to define the new location.
* @param right the vector to add
* @return the new vector that results
*/
point3D & point3D::operator+= (const point3D & right)
{
	for(int d=0; d < TURB3D_DIM ; d++)
		xyz[d] += right.xyz[d] ;
	return *this ;
}

/** Subtract a point to define the difference vector.
* @param right 
* @return the difference vector that results
*/
point3D & point3D::operator-= (const point3D & right)
{
	for(int d=0; d < TURB3D_DIM ; d++)
		xyz[d] -= right.xyz[d] ;
	return *this ;
}

/** Add two 3D points in the sense that the first is a point, the 2nd a vector for translation.
* @param left the point/vector left to the summation sign
* @param right the point/vector right of the summation sign
* @return the sum. This contains the sum of the two inputs in each component.
*/
point3D operator+( const point3D & left, const point3D & right)
{
	return point3D( left.xyz[0]+right.xyz[0], left.xyz[1]+right.xyz[1], left.xyz[2]+right.xyz[2]) ;
}

/** Subtract two 3D points in the sense that the first is a point, the 2nd a point and the result is the vector from the second to the first.
* @param left the point/vector left to the subtraction sign
* @param right the point/vector right of the subtraction sign
* @return The difference.
*/
point3D operator-( const point3D & left, const point3D & right)
{
	return point3D( left.xyz[0]-right.xyz[0], left.xyz[1]-right.xyz[1], left.xyz[2]-right.xyz[2]) ;
}

/** Multiply length by a factor.
* @param c multiplier.
* @return stretched (c >1) or shrinked (c<1)  or reverted (c<0) vector.
*/
point3D operator*( const double c, const point3D & right)
{
	return point3D( c*right.xyz[0], c*right.xyz[1], c*right.xyz[2]) ;
}

/** Print a Position.
* @param os the output stream to print to
* @param some the term to be printed
*/
ostream & operator<<(ostream &os, const point3D & some)
{
	cout << " ( " << some.xyz[0] << " " << some.xyz[1] << " " << some.xyz[2] << " ) " ;
	return os ;
}

/** Ctor 
* @param c the coeffient in front
* @param exx the power of x
* @param exy the power of y
* @param exz the power of z
*/
trino3D::trino3D(const double c, const int exx, const int exy, const int exz) : coef(c)
{
	expo[0] = exx ;
	expo[1] = exy ;
	expo[2] = exz ;
}

/** Evaluation a a specific point in 3D space
*/
double trino3D::at(const point3D & pt) const
{
#ifdef HAVE_GSL
	return coef*gsl_sf_pow_int( pt.xyz[0], expo[0]) *gsl_sf_pow_int( pt.xyz[1], expo[1]) 
			*gsl_sf_pow_int( pt.xyz[2], expo[2]) ;
#else
	return coef*pow( pt.xyz[0], expo[0]) *pow( pt.xyz[1], expo[1]) *pow( pt.xyz[2], expo[2]) ;
#endif
}

/** Gradient evaluation at a specified point in 3D space.
* @param pt the point at which the gradient is computed
* @param gr on return the three components of the gradient
*/
void trino3D::gradat(const point3D & pt, double gr[3]) const
{
#if 0
	for(int i=0 ; i < 3 ; i++)
	{
		gr[i] =  (expo[i] > 0) ? 
			coef*pow( pt.xyz[0], expo[0]) *pow( pt.xyz[1], expo[1]) *pow( pt.xyz[2], expo[2]) ;
			: 0.0 ;
	}
#endif
	gr[0] =  (expo[0] > 0) ? 
			coef*expo[0]*pow( pt.xyz[0], expo[0]-1) *pow( pt.xyz[1], expo[1]) *pow( pt.xyz[2], expo[2])
			: 0.0 ;
	gr[1] =  (expo[1] > 0) ? 
			coef*expo[1]*pow( pt.xyz[0], expo[0]) *pow( pt.xyz[1], expo[1]-1) *pow( pt.xyz[2], expo[2])
			: 0.0 ;
	gr[2] =  (expo[2] > 0) ? 
			coef*expo[2]*pow( pt.xyz[0], expo[0]) *pow( pt.xyz[1], expo[1]) *pow( pt.xyz[2], expo[2]-1)
			: 0.0 ;
}

/** Test whether the polynomial is of some type of given exponents.
* Check whether the current polynomial is compatible with respect to addition.
* @param p the three non-negative exponents of the reference type
* @return true of the current type matches all three exponents of the reference type.
*/
bool trino3D::istype( const int e[3]) const
{
	return (e[0] == expo[0] && e[1] == expo[1] && e[2] == expo[2] );
}

/** Test whether the polynomial is of some type of given exponents.
* Check whether the current polynomial is compatible with respect to addition.
* @param oth the trinomial to match against
* @return true of the current type matches all three exponents of the reference type.
*/
bool trino3D::istype( const trino3D & oth) const
{
	/* forward the test to the more basic implementation with 
	* the function signature of three integer exponents.
	*/
	return istype(oth.expo) ;
}

/** Multiply by a constant
* @param c the constant to multply with
* @return the product that results
*/
trino3D & trino3D::operator*= (const double c)
{
	coef *= c ;
	return *this ;
}

/** Multiply two 3D trinomials.
* @param left the trinomial left to the multiplication sign
* @param right the trinomial right of the multiplication sign
* @return the product. The coefficient is the product of the coeffienst and
*        the exponents are the sums of the individual exponents of x, y and z.
*/
trino3D operator*( const trino3D & left, const trino3D & right)
{
	return trino3D( left.coef*right.coef, left.expo[0]+right.expo[0], left.expo[1]+right.expo[1], left.expo[2]+right.expo[2]) ;
}

/** Print a trinomial term
* @param os the output stream to print to
* @param some the term to be printed
*/
ostream & operator<<(ostream &os, const trino3D & some)
{
	cout << some.coef << " [" << some.expo[0] << "," << some.expo[1] << "," << some.expo[2] << "] " ;
	return os ;
}


Zern3DRadi::Zern3DRadi(const int lowIdx, const int upIdx) : n(lowIdx), l(upIdx)
{
	alpha=(n-l)/2 ;
	double gfac = sqrt(2.*n+3.)/gsl_sf_pow_int(2.,n-l)/gsl_sf_choose((unsigned)n,(unsigned)l) ;
	if ( GSL_IS_ODD(alpha)) 
		gfac *= -1. ;
	for(int s=0 ; s <= alpha; s++)
	{
		double c = gsl_sf_choose((unsigned)n,(unsigned)(alpha-s))
			* gsl_sf_choose((unsigned)(l+s),(unsigned)l)
			* gsl_sf_choose((unsigned)(l+1+n+2*s),(unsigned)(n-l)) ;
		if ( GSL_IS_ODD(s) )
			c *= -1. ;
		coef.push_back(gfac*c) ;
	}
}

#if INCLUDE_SUPERFLUOUS
/** Evaluate at some point 0<=r<=1
* @param r the distance to the origin
* @return R_n^l(r)
* @warn This is not needed anywhere in this program.
*/
double Zern3DRadi::at(double r)
{
	double res=0. ;
	double rpow = gsl_sf_pow_int(r,l) ;
	for(int s=0 ; s <= alpha; s++)
	{
		res += coef[s]*rpow ;
		rpow *= r*r ;
	}
	return res ;
}
#endif


/** Default ctor.
*/
poly3D::poly3D()
{
	/* nothing to do. Empty list of elements */
}

/** Constructor of a radial polynomial r^nl with r^2=x^2+y^2+z^2
* @param nl the power of r. Must be an even integer.
* @warn nl must be an even integer which is not checked in this revision.
*/
poly3D::poly3D(int nl)
{
	/* need only half of this for the implementation */
	nl /= 2 ;

	/* the factorial (n-l)/2 ! 
	*/
	const double nlfac = gsl_sf_fact(nl) ;
	for(int sig1=0; sig1 <= nl ; sig1 ++)
	for(int sig2=0; sig2 <= nl-sig1 ; sig2 ++)
	{
		/* the coefficient is (n-l)/2 !/[sig1 ! sig2 ! (nl/2-sig1-sig2)!]
		*/
		const double c = gsl_sf_fact(sig1)*gsl_sf_fact(sig2)*gsl_sf_fact(nl-sig1-sig2) ;
		const trino3D tmp(nlfac/c, 2*sig1, 2*sig2, 2*(nl-sig1-sig2)) ;

		/* This trinomial is stored into the vector of terms. There is
		* no need tocheck for duplicates because the summation above
		* cannot produce some... the exponents are unique for each pair of 
		* indices in the double sum.
		*/
#if 0
		cout << __LINE__ << " " << 2*nl << " " << tmp << endl ;
#endif
		compo.push_back(tmp) ;
	}
}

/** Constructor for a vector Harmonics r^l Y_l^m.
* @param l the angular momentum quantum number
* @param m the magnetic quantum number in the range 0<=m<=l
* @param cnots true if this is a cosine type, false if a sine type
* @warn there is no check that m=0 is only used with cnots equal to true
*/
poly3D::poly3D(const int l, const int m, const bool cnots)
{
	/* global factor, including the sqrt(2) for the definition
	* of the real-valued Ylm if m<>0
	*/
	double gfac = sqrt((l+0.5)/M_PI*gsl_sf_fact(l-m)*gsl_sf_fact(l+m)) ;

	/* multiply by (-1/2)^m, and recover the sqrt(2) factor if m=0
	*/
	if ( m )
		gfac *= gsl_sf_pow_int(-0.5,m) ;
	else
		gfac *= M_SQRT1_2 ;

	const int lmhalf = (l-m)/2 ;
	for(int sig1=0; sig1 <= lmhalf ; sig1 ++)
	for(int sig2=0; sig2 <= lmhalf-sig1 ; sig2 ++)
	{
		/* the other factor is (-1/4)^(sig1+sig2)/sig1!/sig2!/(m+sig1+sig2)!/(l-m-2sig1-2sig2)!
		*/
		const double fac2 = gsl_sf_pow_int(-0.25,sig1+sig2)/
			(gsl_sf_fact(sig1)*gsl_sf_fact(sig2)*gsl_sf_fact(m+sig1+sig2)*gsl_sf_fact( l-m-2*(sig1+sig2) ));

		/* loop over even j if cosine, over odd j if sine
		*/
		int jstrt = ( cnots ) ? 0 : 1 ;
		for (int j= jstrt ; j <= m ; j += 2 )
		{
			double c = gfac*fac2* gsl_sf_choose((unsigned)m,(unsigned)j) ;
			if ( GSL_IS_ODD(j/2) ) 
				c *= -1 ;
			const trino3D tmp(c, m-j+2*sig1, j+2*sig2, l-m-2*(sig1+sig2)) ;

			/* This trinomial is stored into the vector of terms. There is
			* no need tocheck for duplicates because the summation above
			* cannot produce some... the exponents are unique for each pair of 
			* indices in the double sum.
			*/
#if 0
			cout << __LINE__ << " " << 2*lmhalf << " " << tmp << endl ;
#endif
			compo.push_back(tmp) ;
		}
	}
}

/** Multiply by a constant
* @param c the constant to multply with
* @return the product that results
*/
poly3D & poly3D::operator*= (const double cof)
{
	if ( cof == 0.)
		compo.clear() ;
	else
		for(int c=0 ; c < compo.size() ; c++)
			compo[c] *= cof ;
	return *this ;
}

/** Constructor for the radial Zernike polynomials R_n^l(r)/r^l.
* @param n the main power
* @param l the angular momentum quantum number
*/
poly3D::poly3D(const int n, const int l)
{
	/* construct the R_n^l radial polynomial in front of Y_lm
	*/
	Zern3DRadi R(n,l) ;

	/* accumulate the powers r^(2s) into trinomials
	*/
	for(int s=0 ; s <= R.alpha; s++)
	{
		poly3D r2s(2*s) ;
		r2s *= R.coef[s] ;
		*this += r2s ;
	}
#if 0
	cout << __LINE__ << " " ;
	for(int i=0; i < compo.size() ;i++)
	cout << compo[i] ;
	cout << endl ;
#endif
}

/** Constructor for a general Zernike term R_n^l Y_lm
* @param n the main power
* @param l the angular momentum quantum number
* @param m the magnetic quantum number in the range 0<=m<=l
* @param cnots true if this is a cosine type, false if a sine type
* @warn there is no check that m=0 is only used with cnots equal to true
*/
poly3D::poly3D(const int n, const int l, const int m, const bool cnots)
{
	/* first construct the vector Harmonics r^l Y_lm
	*/
	*this = poly3D(l,m,cnots) ;

	/* construct R_n^l without the r^l factor and multiply with the existing object
	*/
	*this *= poly3D(n,l) ;
}

/** Evaluate it at a point in 3D
* param pt the location of the 3D point
* return the value, which is the sum over all terms.
*/
double poly3D::at(const point3D & pt) const
{
	/* If the current polynomial has no terms, its value is zero.
	*/
	double res = 0. ;
	for(int c =0; c < compo.size() ; c++)
		res += compo[c].at(pt) ;
	return res ;
}

/** Gradient evaluation at a specified point in 3D space.
* @param pt the point at which the gradient is computed
* @param gr on return the three components of the gradient
*/
void poly3D::gradat(const point3D & pt, double gr[3]) const
{
	/* If the current polynomial has no terms, its value is zero.
	*/
	for(int d=0 ; d < TURB3D_DIM ; d++)
		gr[d] =0. ;

	for(int c =0; c < compo.size() ; c++)
	{
		double tgrad[3] ;
		compo[c].gradat(pt,tgrad) ;
		for(int d=0 ; d < TURB3D_DIM ; d++)
			gr[d] +=  tgrad[d] ;
	}
}

/** Test whether a trinomial of an exponential signature is
* already one of the terms.
* @param e the three nonzero exponents of the reference
* @return the index of the first term of that signature, or -1 if not found.
*/
int poly3D::hastype(const int e[TURB3D_DIM]) const
{
	for(int c=0; c < compo.size() ; c++)
		if ( compo[c].istype(e) )
			return c ;
	return -1 ;
}

/** Test whether a trinomial of an exponential signature is
* already one of the terms.
* @param t the trinomial of the reference
* @return the index of the first term of that signature, or -1 if not found.
*/
int poly3D::hastype(const trino3D & t) const
{
	return hastype(t.expo) ;
}

/** Add another trinomial term to this one
* @param oth the trinomial on the right hand side of the equation.
* @todo allow for annihilation of coefficients and elimination of associated products..
*/
poly3D & poly3D::operator += (const trino3D & oth)
{
	/* Check whether the exponent's signature is already present 
	*/
	const int idx = hastype( oth ) ;
	if ( idx >=0 )
		compo[idx].coef += oth.coef ;
	else
		compo.push_back( oth ) ;
	return *this ;
}

/** Multiply all components with the other trinomial term.
* @param oth the trinomial on the right hand side of the equation.
*/
poly3D & poly3D::operator *= (const trino3D & t)
{
	/* If the other term has a coefficient of zero, erase the current one, so
	* it effectively evaluates to zero.
	*/
	if ( t.coef == 0.) 
		compo.clear() ;
	else
	{
		for(int c=0 ; c < compo.size() ; c++)
		{
			/* multiply coefficients */
			compo[c].coef *= t.coef ;
			/* add exponents */
			compo[c].expo[0]  += t.expo[0]  ;
			compo[c].expo[1]  += t.expo[1]  ;
			compo[c].expo[2]  += t.expo[2]  ;
		}
	}
	return *this ;
}

/** Add another polynomial to this one
* @param oth the polynomial on the right hand side of the equation.
*/
poly3D & poly3D::operator += (const poly3D & oth)
{
	/* Add the other polynomial's terms one by one
	*/
	for(int c=0 ; c < oth.compo.size() ; c++)
		*this += oth.compo[c] ;
	return *this ;
}

/** Multiply by another polynomial.
* @param oth the polynomial on the right hand side of the equation.
*/
poly3D & poly3D::operator *= (const poly3D & oth)
{
	/* No terms in the other polynomial means the result is zero
	*/
	if ( oth.compo.size() ) 
	{
		/* Start with a polynomial of value zero, no terms
		*/
		poly3D res ;
	
		/* Build with the distributive law all pairs of products and accumulate them
		*/
		for (int i=0 ; i < compo.size() ; i++)
		for (int j=0 ; j < oth.compo.size() ; j++)
			res += compo[i]*oth.compo[j] ;

		/* overwrite the present polynomial */
		compo = res.compo ;
	}
	else
		compo.clear() ;
	return *this ;
}

/** Multiply two 3D polynomials.
* @param left the polynomial left to the multiplication sign
* @param right the polynomial right of the multiplication sign
* @return the polynomial which is the sum over all products of the components.
*/
poly3D operator*( const poly3D & left, const poly3D & right)
{
#if 0
	poly3D res(left) ;
	return (res *= right );
#endif

	/* Start with a polynomial of value zero, no terms
	*/
	poly3D res ;

	/* Build with the distributive law all pairs of products and accumulate them
	*/
	for (int i=0 ; i < left.compo.size() ; i++)
	for (int j=0 ; j < right.compo.size() ; j++)
		res += left.compo[i]*right.compo[j] ;
	return res ;
}

/**
*/
ostream & operator<<(ostream &os, const poly3D & some)
{
	for(int i=0; i < some.compo.size() ;i++)
		cout << some.compo[i] ;
	cout << endl ;
	return os ;
}
