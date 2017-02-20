/* jacobi.c
 * Copyright (C) 2006 Paulo José Saiz Jabardo
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

/* Author:  Paulo J. Saiz Jabardo */






/** \file jacobi.c
 
 \brief Calculation of Jacobi Polynomials

 This file implements the function that calculates Jacobi polynomials and its derivatives.
 It tries to follow the interfaces used in the GNU Scientific library. The estimated errors
 Have not been implemented yet.

  The main sources used to develop this library where:

  (1) Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
      Dover Publications 1965.

  (2) Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
      Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.


   The Jacobi polynomials of order n are defined as follows

   \f[ P_n^{\alpha,\beta}(x) \f]

   The calculation uses the recursive definition of the Jacobi polynomial:

   \f[
   a^1_n P_{n+1}^{\alpha,\beta}(x) = (a^2_n + a^3_n x)P_n^{\alpha,\beta}(x) - a^4_n P_{n-1}^{\alpha,\beta}(x)
   \f]
   
   In this function, \f$-1\le x \le 1\f$ is the point where the function is to be calculated, n is the order of the
   polynomial and \f$a^i_n\f$ are coefficients given by

   \f[
     a^1_n = 2(n+1)(n+\alpha+\beta+1)(2n+\alpha+\beta)
   \f]
   \f[
     a^2_n = (2n+\alpha+\beta+1)(\alpha^2-\beta^2)
   \f]
   \f[
     a^3_n = (2n + \alpha + \beta)(2n + \alpha + \beta + 1)(2n + \alpha +\beta + 2)
   \f]
   \f[
     a^4_n = 2(n+\alpha)(n+\beta)(2n+\alpha+\beta+2)
   \f]



   There is also a recursive definition of the derivative of Jacobi polynomials but in this library, the derivatives
   of Jacobi polynomials will be defined in terms of Jacobi polynomials:

   \f[
     \frac{d}{dx}P^{\alpha,\beta}_n(x) = \frac{1}{2}(\alpha + \beta + n + 1) P_n^{\alpha+1,\beta+1}(x)
   \f]

   Chebychev polynomials are a special case of Jacobi polynomials with \f$\alpha,\beta=-1/2\f$:

   \f[
      T_n(x) = \frac{2^{2n}(n!)^2}{(2n)!} P_n^{-1/2,-1/2}(x)
   \f]

   Since the Zeros of Chebychev are explicitly known, an iteration using these zeros as initial guess is used
   to calculate zeros of generic Jacobi polynomials.
   

 */    
 
 


#include <stdlib.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>



#include "jacobi.h"





  

/** Calculates the 0th order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param  result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
  
 */
int
jac_jacobi_P0_e (double x, double a, double b, gsl_sf_result *result)
{
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
}

/** Calculates the 1st order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
  
 */
int
jac_jacobi_P1_e (double x, double a, double b, gsl_sf_result *result)
{
    result->val = 0.5 * (a - b + (a + b + 2.0) * x);
    result->err = 0.0;
    return GSL_SUCCESS;
}



/** Calculates the nth order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param n An integer specifying the order of the polynomial
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
 
 */
int
jac_jacobi_e (double x, int n, double a, double b, gsl_sf_result *result)
{
    
    if (n==0)
    {
	result->val = 1.0;
	result->err = 0.0;
	return GSL_SUCCESS;
    }
    else if (n==1)
    {
	result->val = 0.5 * (a - b + (a + b + 2.0)*x);
	result->err = 0.0;
	return GSL_SUCCESS;
    }
    else
    {
	
	double p0, p1, a1, a2, a3, a4, p2=0.0;
	int i;
	p0 = 1.0;
	p1 = 0.5 * (a - b + (a + b + 2)*x);
    
	for(i=1; i<n; ++i){
	    a1 = 2.0*(i+1.0)*(i+a+b+1.0)*(2.0*i+a+b);
	    a2 = (2.0*i+a+b+1.0)*(a*a-b*b);
	    a3 = (2.0*i+a+b)*(2.0*i+a+b+1.0)*(2.0*i+a+b+2.0);
	    a4 = 2.0*(i+a)*(i+b)*(2.0*i+a+b+2.0);
	    p2 = 1.0/a1*( (a2 + a3*x)*p1 - a4*p0);
	
	    p0 = p1;
	    p1 = p2;
	}
    
	result->val = p2;
	result->err = 0.0;
	return GSL_SUCCESS;
    }
}


/** Explicit calculation of the 0th order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_jacobi_P0(double x, double a, double b)
{
    return 1.0;
}

/** Explicit calculation of the 1st order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_jacobi_P1(double x, double a, double b)
{
    return 0.5 * (a - b + (a + b + 2.0) * x);
}




/** Explicit calculation of the nth order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param n An integer specifying the order of the polynomial
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_jacobi (double x, int n, double a, double b)
{
    
    if (n==0)
    {
	return 1.0;
    }
    else if (n==1)
    {
	return  0.5 * (a - b + (a + b + 2.0)*x);
    }
    else
    {
	
	double p0, p1, a1, a2, a3, a4, p2=0.0;
	int i;
	p0 = 1.0;
	p1 = 0.5 * (a - b + (a + b + 2)*x);
    
	for(i=1; i<n; ++i){
	    a1 = 2.0*(i+1.0)*(i+a+b+1.0)*(2.0*i+a+b);
	    a2 = (2.0*i+a+b+1.0)*(a*a-b*b);
	    a3 = (2.0*i+a+b)*(2.0*i+a+b+1.0)*(2.0*i+a+b+2.0);
	    a4 = 2.0*(i+a)*(i+b)*(2.0*i+a+b+2.0);
	    p2 = 1.0/a1*( (a2 + a3*x)*p1 - a4*p0);
	
	    p0 = p1;
	    p1 = p2;
	}
    
	return p2;
    }
}





/** Calculates the derivative of the 0th order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
  
 */
int
jac_djacobi_P0_e (double x, double a, double b, gsl_sf_result *result)
{
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
}

/** Calculates the derivative of the 1st order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
  
 */
int
jac_djacobi_P1_e (double x, double a, double b, gsl_sf_result *result)
{
    result->val = 0.5 * (a + b + 2.0);
    result->err = 0.0;
    return GSL_SUCCESS;
}


/** Calculates the derivative of the nth order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param n Integer specifying the order of the polynomial
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \param result A pointer to a gsl_sf_result structure where the result and estimated error is stored
  \return GSL_SUCCESS if successful
  
 */
int
jac_djacobi_e (double x, int n, double a, double b, gsl_sf_result *result)
{
    int code = jac_jacobi_e(x, n-1, a+1.0, b+1.0, result);
    
    result->val *= 0.5  * (a + b + n + 1.0);
    result->err = 0.0;
    return code;

}


/** Explicit calculation of the derivative of the 0th order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_djacobi_P0(double x, double a, double b)
{
    return 0.0;
}

/** Explicit calculation of the derivative of the 1st order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_djacobi_P1(double x, double a, double b)
{
    return 0.5 * (a + b + 2.0);
}


/** Explicit calculation of the derivative of the nth order Jacobi polynomial

  \param x The point where the function should be calculated. \f$-1\le x \le 1\f$
  \param n The order of the polynomial
  \param a \f$\alpha\f$ parameter of the Jacobi Polynomial
  \param b \f$\beta\f$ parameter of the Jacobi Polynomial
  \return The calculated value of the Jacobi polynomial
  
 */
double
jac_djacobi (double x, int n, double a, double b)
{
    
    if (n==0)
    {
	return 0.0;
    }
    else if (n==1)
    {
	return  0.5 * (a + b + 2.0);
    }
    else
    {
	return 0.5 * (a + b + n + 1.0) * jac_jacobi(x, n-1, a+1.0, b+1.0);
    }
}



/** Calculates the nth order Jacobi polynomials at an array of np points given by the vector x.
    This function needs two workspaces of size np to do its calculation. If the first one (ws1) is
    NULL, the function will allocate the necessary memory (using malloc) to do the calculation

    \param np Number of points where the polynomial will be calculated
    \param x A vector of length np containing the points where the polynomials should be calculated. \f$-1\le x_i \le 1\f$.
    \param n Order of Jacobi polynomial to be calculated
    \param result_array An array of length np where the result will be stored
    \param a \f$\alpha\f$ parameter of Jacobi polynomial
    \param b \f$\beta\f$ parameter of Jacobi polynomial
    \param ws A pointer to a block of memory at least 2*np doubles long to be used as workspace. If it is NULL, memory will be allocated using malloc and released at the end
    \return Error code defined in gsl_errno.h or GSL_SUCCESS if everything was fine
 */
int 
jac_jacobi_array (int np, const double *x, int n, double * result_array,
		  double a, double b, double *ws)
{



    int i;
    // Computes the jacobi polyniomials on several points
    if (n == 0)
    {
	for(i = 0; i < np; ++i)
	    result_array[i] = 1.0;
	return GSL_SUCCESS;
    }
    
    if (n == 1)
    {
	for (i = 0; i < np; ++i)
	    result_array[i] = 0.5 * (a - b + (a + b + 2.0)*x[i]);
	return GSL_SUCCESS;
    }
    
    // General case,
    double *pnm1, *pnm2;
    int mem_allocated=0;
    if (ws==NULL)
    {
	pnm1 = (double *) malloc(2*np*sizeof(double));
	if (!pnm1) return GSL_ENOMEM;
	
	mem_allocated = 1;
	pnm2 = pnm1 + np;
    } else
    {
	pnm1 = ws;
	pnm2 = ws + np;
    }
    
    for(i = 0; i < np; ++i){
	pnm2[i] = 1.0;
	pnm1[i] = 0.5 * (a - b + (a + b + 2.0)*x[i]);
    }
    
    // Start iterating:
    int k;
    double a1, a2, a3, a4;
    for(k=1; k<n; ++k){
	
	a1 = 2.0 * (k + 1.0) * (k + a + b + 1.0) * (2.0 * k + a + b);
	a2 = (2.0 * k + a + b + 1.0) * (a * a - b * b) / a1;
	a3 = (2.0 * k + a + b) * (2.0 * k + a + b + 1.0) * (2.0 * k + a + b + 2.0) / a1;
	a4 = 2.0 * (k + a) * (k + b) * (2.0 * k + a + b + 2.0) / a1;
	for (i = 0; i < np; ++i)
	{
	    result_array[i] = (a2 + a3 * x[i]) * pnm1[i] - a4 * pnm2[i];
	    pnm2[i] = pnm1[i];
	    pnm1[i] = result_array[i];
	}
    }


    if (mem_allocated)
    {
	free(pnm1);
    }
    
    return GSL_SUCCESS;
    

}








/** Calculates the derivative of the nth order Jacobi polynomial at an array of np points given by the vector x.
    This function needs two workspaces of size np to do its calculation. If the first one (ws1) is
    NULL, the function will allocate the necessary memory (using malloc) to do the calculation

    \param np Number of points where the polynomial will be calculated
    \param x A vector of length np containing the points where the polynomials should be calculated. \f$-1\le x_i \le 1\f$.
    \param n Order of the Jacobi polynomial whose derivative should be calculated
    \param result_array An array of length np where the result will be stored
    \param a \f$\alpha\f$ parameter of Jacobi polynomial
    \param b \f$\beta\f$ parameter of Jacobi polynomial
    \param ws A pointer to a block of memory 2*np doubles long to be used as workspace. If it is NULL, memory will be allocated using malloc and released at the end
    \return Error code defined in gsl_errno.h or GSL_SUCCESS if everything was fine
 */
int
jac_djacobi_array(int np, const double *x, int n, double * result_array, double a, double b, double *ws)
{
    int code = jac_jacobi_array(np, x, n-1, result_array, a+1.0, b+1.0, ws);
    if (code)
	return code;
    
    int i;
    for (i = 0; i < np; ++i)
	result_array[i] *= 0.5*(a+b+n+1.0);
    return code;
}






/** This function uses an iteration technique to calculate the zeros of a Jacobi polynomial:

  \f[
  P_m^{\alpha,\beta}(x_i) = 0
  \f]

  The algorithm was adapted from the book by Karniadakis and Sherwin. It uses as an initial
  guess the zeros of Chebyshev polynomials.

  \param x Pointer to an array where the zeros of the Jacobi polynomial should be stored
  \param m Order of the polynomial
  \param alpha \f$\alpha\f$ parameter of the Jacobi polynomial
  \param beta \f$\beta\f$ parameter of the Jacobi polynomial

 */
int
jac_jacobi_zeros(double *x, int m, double alpha, double beta)
{
    // This function returns all of the zeros of a jacobi polynomial. Algorithm taken from B.2 [1];
    
    // A good estimate for the roots is the zero for the chebyshev
    //polynomial which is readily known
    
    int k, i;
    double r, delta, s, poly;
    const double EPS = 100*GSL_DBL_EPSILON;
    const int MAXITER = 200;
    int iter;
    for (k = 0; k < m; k++)
    {
	// Initial guess
	r = -cos( (2*k+1.0)/(2*m) * M_PI);
	if (k > 0) r = (r + x[k-1]) / 2;
	iter = 0;
	do
	{
	    
	    s = 0;
	    
	    for (i=0; i < k; i++)
		s += 1 / (r - x[i]);
	    
	    poly = jac_jacobi(r, m, alpha, beta);
	    delta = - poly / (jac_djacobi(r, m, alpha, beta) - poly*s);
	    r += delta;
	    ++iter;
	    if (iter > MAXITER)
		return GSL_CONTINUE;
	    
	} while (fabs(delta) > EPS);
	
	x[k] = r;
    }
    return GSL_SUCCESS; 
}


