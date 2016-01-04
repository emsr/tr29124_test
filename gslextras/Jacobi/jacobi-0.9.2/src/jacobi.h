/* jacobi.h
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






/** \file jacobi.h
  \brief Header file for jacobi library
  
  The jacobi library implements Jacobi polynomials and Gauss-Jacobi quadrature related functions.
  This library was developed to be used in high order finite element methods (Spectral Element/hp
  methods. The interfaces were inspired by the software Polylib (http://www.nektar.info).

  The main sources used to develop this library where:

  (1) Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
      Dover Publications 1965.

  (2) Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
      Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.

  This software is supposed to be an add-on module for the GNU Scientific Library
*/






#ifndef __JACOBI_H__
#define __JACOBI_H__

#ifdef __cplusplus
//extern "C"{
#endif

#include <gsl/gsl_sf_result.h>



/** \brief Enumeration with differing types of quadrature defined in the library

  The jac_quad_type is used to determine the type of quadrature that is being used.
   
 */
enum jac_quad_type{
  JAC_GJ,   ///< Gauss-Jacobi quadrature
  JAC_GLJ,  ///< Gauss-Lobatto-Jacobi quadrature
  JAC_GRJM, ///< Gauss-Radau-Jacobi quadrature including the node -1
  JAC_GRJP  ///< Gauss-Radau-Jacobi quadrature including the node +1
};


/** \brief Struct used to store information and memory on the quadrature
    
  This strucuture is used to store data about the quadrature. It can be later used to integrate, 
  derive or interpolate functions. It should be created using the jac_quadrature_alloc function.

*/
typedef struct{ 
  int Q;   			///< Number of quadrature points
  enum jac_quad_type type;      ///< Quadrature type
  double alpha;                 ///< Alpha weight of the quadrature
  double beta;                  ///< Beta weight of the quadrature
  
  double *x;                    ///< Array that stores the nodes coordinates
  double *w;                    ///< Array that stores quadrature weights
  double *D;                    ///< Array that stores the derivative matrix
  double *Imat;                 ///< Array That stores the interpolation matrix
  int np;                       ///< Number of interpolation points
  double *xp;                   ///< Interpolation points
  
} jac_quadrature;



		 
		 

/// Calculates the Jacobi polynomial of order 0
int jac_jacobi_P0_e (double x, double a, double b, gsl_sf_result *result);
/// Calculates the Jacobi polynomial of order 1
int jac_jacobi_P1_e (double x, double a, double b, gsl_sf_result *result);
/// Calculates the Jacobi polynomial of order n
int jac_jacobi_e (double x, int n, double a, double b, gsl_sf_result *result);
/// Calculates the Jacobi polynomial of order 0
double jac_jacobi_P0(double x, double a, double b);
/// Calculates the Jacobi polynomial of order 1
double jac_jacobi_P1(double x, double a, double b);
/// Calculates the Jacobi polynomial of order n
double jac_jacobi (double x, int n, double a, double b);
/// Calculates the derivative of the Jacobi polynomial of order 0
int jac_djacobi_P0_e (double x, double a, double b, gsl_sf_result *result);
/// Calculates the derivative of the Jacobi polynomial of order 1
int jac_djacobi_P1_e (double x, double a, double b, gsl_sf_result *result);
/// Calculates the derivative of the Jacobi polynomial of order n
int jac_djacobi_e (double x, int n, double a, double b, gsl_sf_result *result);
/// Calculates the derivative of the Jacobi polynomial of order 0
double jac_djacobi_P0(double x, double a, double b);
/// Calculates the derivative of the Jacobi polynomial of order 1
double jac_djacobi_P1(double x, double a, double b);
/// Calculates the derivative of the Jacobi polynomial of order n
double jac_djacobi (double x, int n, double a, double b);
/// Calculates Jacobi polynomials at an array of points
int jac_jacobi_array (int np, const double *x, int n, double * result_array,
		      double a, double b, double *ws);
/// Calculates the derivative of Jacobi polynomials at an array of points
int jac_djacobi_array(int np, const double *x, int n, double * result_array,
		      double a, double b, double *ws);
/// Calculates the zeros of Jacobi polynomials in the interval -1 up to 1
int jac_jacobi_zeros(double *x, int m, double alfa, double beta);


/// Calculates the integral of a function with values at quadrature points given by f with quadrature weights w
double jac_integrate(jac_quadrature *quad, double *f);
/// Calculates the derivative  at quadrature points of a function f and derivative matrix D
int jac_differentiate(jac_quadrature *quad, double *f, double *d);
/// Interpolates the function given by f using the interpolation matrix
int jac_interpolate(jac_quadrature *quad, double *f, double *fout);

/// Calculates the Gauss-Jacobi quadrature points
int jac_zeros_gj(double *z, const int Q, double alpha, double beta);
/// Calculates the Gauss-Jacobi quadrature weights
int jac_weights_gj(double *z, double *w, const int Q, const double alpha,
		   const double beta, double *ws);
/// Calculates the Gauss-Jacobi Derivative matrix
int jac_diffmat_gj(double *z, double *D, const int Q, double alpha, double beta, double *ws);
/// Calculates the Lagrange polynomials through Gauss-Jacobi Quadratures
double jac_lagrange_gj(int i, double zz, int Q, double *z, double alpha, double beta);

/// Calculates the Gauss-Lobatto-Jacobi quadrature points
int jac_zeros_glj(double *z, const int Q, double alpha, double beta);
/// Calculates the Gauss-Lobatto-Jacobi quadrature weights
int jac_weights_glj(double *z, double *w, const int Q, double alpha, double beta, double *ws1);
/// Calculates the Gauss-Lobatto-Jacobi Derivative matrix
int jac_diffmat_glj(double *z, double *D, const int Q, double alpha, double beta, double *ws);
/// Calculates the Lagrange polynomials through Gauss-Lobatto-Jacobi Quadratures
double jac_lagrange_glj(int i, double zz, int Q, double *z, double alpha, double beta);

/// Calculates the Gauss-Radau-Jacobi quadrature points (point -1 included)
int jac_zeros_grjm(double *z, const int Q, double alpha, double beta);
/// Calculates the Gauss-Radau-Jacobi quadrature weights (point -1 included)
int jac_weights_grjm(double *z, double *w, const int Q, double alpha, double beta, double *ws);
/// Calculates the Gauss-Radau-Jacobi Derivative matrix (point -1 included)
int jac_diffmat_grjm(double *z, double *D, const int Q, double alpha, double beta, double *ws);
/// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi Quadratures (point -1 included)
double jac_lagrange_grjm(int i, double zz, int Q, double *z, double alpha, double beta);

/// Calculates the Gauss-Radau-Jacobi quadrature points (point +1 included)
int jac_zeros_grjp(double *z, const int Q, double alpha, double beta);
/// Calculates the Gauss-Radau-Jacobi quadrature weights (point +1 included)
int jac_weights_grjp(double *z, double *w, const int Q, double alpha, double beta, double *ws);
/// Calculates the Gauss-Radau-Jacobi Derivative matrix (point +1 included)
int jac_diffmat_grjp(double *z, double *D, const int Q, double alpha, double beta, double *ws);
/// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi Quadratures (point +1 included)
double jac_lagrange_grjp(int i, double zz, int Q, double *z, double alpha, double beta);


/// Calculates the interpolation matrix for GJ quadrature
int jac_interpmat_gj(double *imat, double *zp, int np, double *z, int Q, double alpha, double beta);
/// Calculates the interpolation matrix for GLJ quadrature
int jac_interpmat_glj(double *imat, double *zp, int np, double *z, int Q, double alpha, double beta);
/// Calculates the interpolation matrix for GRJM quadrature
int jac_interpmat_grjm(double *imat, double *zp, int np, double *z, int Q, double alpha, double beta);
/// Calculates the interpolation matrix for GRJP quadrature
int jac_interpmat_grjp(double *imat, double *zp, int np, double *z, int Q, double alpha, double beta);



/// Allocates the jac_quadrature data structure to be used later
jac_quadrature *jac_quadrature_alloc(int nq);
/// Calculates the quadrature zeros, weights and derivative matrix
int jac_quadrature_zwd(jac_quadrature *quad, enum jac_quad_type qtype, double a, double b, double *ws);
/// Reclaims memory allocated by jac_quadrature_alloc
void jac_quadrature_free(jac_quadrature *quad);

/// Allocates memory to an interpolation matrix
int jac_interpmat_alloc(jac_quadrature *quad, int npoints, double *xp);
/// Reclaims memory allocated by the jac_interpmat_alloc
void jac_interpmat_free(jac_quadrature *quad);








#ifdef __cplusplus
  //}
#endif

#endif /* __JACOBI_H__ */
