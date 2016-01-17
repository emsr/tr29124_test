#ifndef POLY3D_H
#define POLY3D_H

#include <cmath>
#include <vector>
#include <cstring>
#include <cstdarg>

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std ;

#define TURB3D_DIM 3

/** A point in a 3D domain.
* @since 2008-09-22
*/
class point3D {
public:
	/** the 3 Cartesian coordinates.
	*/
	double xyz[3] ;

#ifdef INCLUDE_SUPERFLUOUS
	/** Sphere radius
	*/
	double rad ;

	point3D( const double cart[TURB3D_DIM], double radius =1.) ;
#endif

	point3D( ) ;

	point3D( const double X, const double Y, const double Z) ;

	point3D( const double cart[TURB3D_DIM]) ;

	point3D( const point3D & orig) ;

	double dist() const ;
	double dist( const point3D &oth) const ;

	void scale(const double radius) ;

	void normalize(const double newleng) ;

	void phithet( double rtp[3]) const ;

	point3D & operator= (const point3D  & right) ;

	point3D & operator*= (const double c) ;

	point3D & operator+= (const point3D & right) ;
	point3D & operator-= (const point3D & right) ;

protected:
private:
} ; /* point3D */

point3D operator+( const point3D & left, const point3D & right) ;
point3D operator-( const point3D & left, const point3D & right) ;

point3D operator*( const double c, const point3D & right) ;

ostream & operator<<(ostream &os, const point3D & some) ;

/** A term in a 3D polynomial
* @since 2008-09-24
*/
class trino3D {
public:
	/** The expansion coefficient
	*/
	double coef ;

	/** The non-negative powers of x, y and z
	*/
	int expo[TURB3D_DIM] ;

	trino3D(const double c, const int exx, const int exy, const int exz) ;

	double at(const point3D & pt) const ;

	void gradat(const point3D & pt, double gr[3]) const ;

	bool istype( const int e[3]) const ;

	bool istype( const trino3D & oth) const ;

	trino3D & operator*= (const double c) ;

protected:
private:
} ;

trino3D operator*( const trino3D & left, const trino3D & right) ;

ostream & operator<<(ostream &os, const trino3D & some) ;

/** Zernike 3D radial polynomial
* @since 2008-09-25
*/
class Zern3DRadi {
public:
	/** main quantum number, lower index
	*/
	int n ;

	/** upper index, representing the family
	*/
	int l ;

	/* coefficients for powers r^(l+2s), first value for s=0, last for s=(n-l)/2
	*/
	vector<double> coef ;

	/** auxiliary mixed excess index
	*/
	int alpha ;

	Zern3DRadi(const int lowIdx, const int upIdx) ;

#if INCLUDE_SUPERFLUOUS
	double at(double r) ;
#endif
protected:
private:
}; /* Zern3DRadi */


/** Polynomial in 3D
* @since 2008-09-24
*/
class poly3D {
public:
	/** The individual terms.
	*/
	vector<trino3D> compo ;

	poly3D() ;

	poly3D(int nl) ;

	poly3D(const int l, const int m, const bool cnots) ;

	poly3D & operator*= (const double cof) ;

	poly3D(const int n, const int l) ;

	poly3D(const int n, const int l, const int m, const bool cnots) ;

	double at(const point3D & pt) const ;

	void gradat(const point3D & pt, double gr[3]) const ;

	int hastype(const int e[TURB3D_DIM]) const ;

	int hastype(const trino3D & t) const ;

	poly3D & operator += (const trino3D & oth) ;

	poly3D & operator *= (const trino3D & t) ;

	poly3D & operator += (const poly3D & oth) ;

	poly3D & operator *= (const poly3D & oth) ;

protected:
private:
} ; /* Poly3D */

poly3D operator*( const poly3D & left, const poly3D & right) ;

ostream & operator<<(ostream &os, const poly3D & some) ;

#endif /* POLY3D_H */
