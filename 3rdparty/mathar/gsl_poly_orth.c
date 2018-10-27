/* gsl/poly/gsl_poly_orth.c
 *
 * Copyright (C) 2006 Richard J. Mathar
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_pow_int.h>

/**
* @param[out] dd divided differences for Hermite Polynomial H(n,x).
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[out] xa abscissa values for the dd array
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[in] size One more than the index n of the polynomial; also the
*        length of the arrays dd and xa.
* A call of the function prepares the arrays dd and xa for calls
* with gsl_poly_dd_eval() or gsl_poly_dd_taylor() .
* @since 2006-09-09
* @author Richard J. Mathar
*/
int gsl_poly_dd_init_orthH (double dd[], double xa[], const size_t size)
{
	int l ;
	const int n=size-1 ;
	double fact= gsl_pow_int(2.,n) ;
	/* We use the representation 25.1.10 of the divided differences,
	* evaluating the derivatives of the Hermite polynomial at x=0.
	*/
	for(l=0 ; l < size ; l++)
		xa[l]=0. ;

	/* Since the polynomial has even parity for even n and odd parity for
	* odd n, the l'th derivative at x=0 is zero for odd n-l.
	*/
	for(l= n-1; l >= 0 ; l -= 2)
		dd[l]=0. ;

	/* the l'th derivative is generally computed by eq 22.3.10 of
	* Abramowitz & Stegun
	*/
	for(l= n; l >= 0 ; l -= 2)
	{
		dd[l]= fact ;
		fact *= -0.5*l*(l-1)/(double)(n-l+2) ;
	}

  	return GSL_SUCCESS;
}

/**
* @param[out] dd divided differences for Laguerre Polynomial L(n,x).
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[out] xa abscissa values for the dd array
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[in] size One more than the index n of the polynomial; also the
*        length of the arrays dd and xa.
* After a call of the function, the arrays dd and xa can be used for calls
* with gsl_poly_dd_eval() or gsl_poly_dd_taylor() .
* @since 2006-09-09
* @author Richard J. Mathar
*/
int gsl_poly_dd_init_orthL (double dd[], double xa[], const size_t size)
{
	int l ;
	const int n=size-1 ;
	double fact= 1.0 ;
	/* We use the representation 25.1.10 of the divided differences,
	* evaluating the derivatives of the Laguerre polynomial at x=0.
	*/
	for(l=0 ; l < size ; l++)
		xa[l]=0. ;

	/* the l'th derivative is generally computed by eq 22.3.9 of
	* Abramowitz & Stegun for the case alpha=0.
	*/
	for(l= 0; l <= n ; l ++)
	{
		dd[l]= fact ;
		fact *= (l-n)/gsl_pow_2(1+l) ;
	}

  	return GSL_SUCCESS;
}

/**
* @param[out] dd divided differences for Chebyshev Polynomial of the 1st kind T(n,x).
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[out] xa abscissa values for the dd array
*   Array must have been allocated with n+1 elements on input and is modified
*   on output.
* @param[in] size One more than the index n of the polynomial; also the
*        length of the arrays dd and xa.
* After a call of the function, the arrays dd and xa can be used for calls
* with gsl_poly_dd_eval() or gsl_poly_dd_taylor() .
* @since 2006-09-09
* @author Richard J. Mathar
*/
int gsl_poly_dd_init_orthT (double dd[], double xa[], const size_t size)
{
	int l ;
	const int n=size-1 ;
	double fact= gsl_pow_int(2.,n-1) ;
	/* We use the representation 25.1.10 of the divided differences,
	* evaluating the derivatives of the Chebyshev polynomial at x=0.
	*/
	for(l=0 ; l < size ; l++)
		xa[l]=0. ;

	/* Since the polynomial has even parity for even n and odd parity for
	* odd n, the l'th derivative at x=0 is zero for odd n-l.
	*/
	for(l= n-1; l >= 0 ; l -= 2)
		dd[l]=0. ;

	/* the l'th derivative is generally computed by eq 22.3.6 of
	* Abramowitz & Stegun, here only at x=0.
	*/
	for(l= n; l >= 0 ; l -= 2)
	{
		dd[l]= fact ;
		fact *= l*(1-l)/(double)((n-l+2)*(n+l-2)) ;
	}

  	return GSL_SUCCESS;
}
