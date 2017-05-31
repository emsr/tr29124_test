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

/** 
 * \brief Enumeration with differing types of quadrature defined in the library
 * The jac_quad_type is used to determine the type of quadrature that is being used.
 */
enum jac_quad_type
{
  Gauss,          ///< Gauss quadrature
  Gauss_Lobatto,  ///< Gauss-Lobatto quadrature
  Gauss_Radau_m1, ///< Gauss-Radau quadrature including the node -1
  Gauss_Radau_p1  ///< Gauss-Radau quadrature including the node +1
};


/**
 * \brief Struct used to store information and memory on the quadrature
 *
 * This strucuture is used to store data about the quadrature. It can be later used to integrate, 
 * derive or interpolate functions. It should be created using the jac_quadrature_alloc function.
 */
template<typename _Tp>
  struct jac_quadrature
  {
    /// Number of quadrature points
    int Q;

    /// Quadrature type
    enum jac_quad_type type;

    /// Alpha weight of the quadrature
    _Tp alpha;

    /// Beta weight of the quadrature
    _Tp beta;

  private:

    /// Array that stores the nodes coordinates
    _Tp* x;

    /// Array that stores quadrature weights
    _Tp* w;

    /// Array that stores the derivative matrix
    _Tp* D;

    /// Array That stores the interpolation matrix
    _Tp* Imat;

    /// Number of interpolation points
    int np;

    /// Interpolation points
    _Tp* xp;
  };


/// Allocates the jac_quadrature data structure to be used later
template<typename _Tp>
  jac_quadrature<_Tp>* jac_quadrature_alloc(int nq);

/// Calculates the quadrature zeros, weights and derivative matrix
template<typename _Tp>
  int jac_quadrature_zwd(jac_quadrature<_Tp>* quad, enum jac_quad_type qtype, _Tp a, _Tp b, _Tp* ws);

/// Reclaims memory allocated by jac_quadrature_alloc
template<typename _Tp>
  void jac_quadrature_free(jac_quadrature<_Tp>* quad);

/// Allocates memory to an interpolation matrix
template<typename _Tp>
  int jac_interpmat_alloc(jac_quadrature<_Tp>* quad, int npoints, _Tp* xp);

/// Reclaims memory allocated by the jac_interpmat_alloc
template<typename _Tp>
  void jac_interpmat_free(jac_quadrature<_Tp>* quad);



/// Calculates the derivative of the Jacobi polynomial of order n
template<typename _Tp>
  _Tp jac_djacobi(_Tp x, int n, _Tp a, _Tp b);

/// Calculates Jacobi polynomials at an array of points
template<typename _Tp>
  int jac_jacobi_array(int np, const _Tp* x, int n, _Tp* result_array,
		      _Tp a, _Tp b, _Tp* ws);

/// Calculates the derivative of Jacobi polynomials at an array of points
template<typename _Tp>
  int jac_djacobi_array(int np, const _Tp* x, int n, _Tp* result_array,
		      _Tp a, _Tp b, _Tp* ws);

/// Calculates the zeros of Jacobi polynomials in the interval -1 up to 1
template<typename _Tp>
  int jac_jacobi_zeros(_Tp* x, int m, _Tp alpha, _Tp beta);


/// Calculates the integral of a function with values at quadrature points given by f with quadrature weights w
template<typename _Tp>
  _Tp jac_integrate(jac_quadrature<_Tp>* quad, _Tp* f);

/// Calculates the derivative at quadrature points of a function f and derivative matrix D
template<typename _Tp>
  int jac_differentiate(jac_quadrature<_Tp>* quad, _Tp* f, _Tp* d);

/// Interpolates the function given by f using the interpolation matrix
template<typename _Tp>
  int jac_interpolate(jac_quadrature<_Tp>* quad, _Tp* f, _Tp* fout);


/// Calculates the Gauss-Jacobi quadrature points
template<typename _Tp>
  int jac_zeros_gj(_Tp* z, const int Q, _Tp alpha, _Tp beta);

/// Calculates the Gauss-Jacobi quadrature weights
template<typename _Tp>
  int jac_weights_gj(_Tp* z, _Tp* w, const int Q, const _Tp alpha,
		   const _Tp beta, _Tp* ws);

/// Calculates the Gauss-Jacobi Derivative matrix
template<typename _Tp>
  int jac_diffmat_gj(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Lagrange polynomials through Gauss-Jacobi Quadratures
template<typename _Tp>
  _Tp jac_lagrange_gj(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta);

/// Calculates the interpolation matrix for GJ quadrature
template<typename _Tp>
  int jac_interpmat_gj(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta);


/// Calculates the Gauss-Lobatto-Jacobi quadrature points
template<typename _Tp>
  int jac_zeros_glj(_Tp* z, const int Q, _Tp alpha, _Tp beta);

/// Calculates the Gauss-Lobatto-Jacobi quadrature weights
template<typename _Tp>
  int jac_weights_glj(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws1);

/// Calculates the Gauss-Lobatto-Jacobi Derivative matrix
template<typename _Tp>
  int jac_diffmat_glj(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Lagrange polynomials through Gauss-Lobatto-Jacobi Quadratures
template<typename _Tp>
  _Tp jac_lagrange_glj(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta);

/// Calculates the interpolation matrix for GLJ quadrature
template<typename _Tp>
  int jac_interpmat_glj(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta);


/// Calculates the Gauss-Radau-Jacobi quadrature points (point -1 included)
template<typename _Tp>
  int jac_zeros_grjm(_Tp* z, const int Q, _Tp alpha, _Tp beta);

/// Calculates the Gauss-Radau-Jacobi quadrature weights (point -1 included)
template<typename _Tp>
  int jac_weights_grjm(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Gauss-Radau-Jacobi Derivative matrix (point -1 included)
template<typename _Tp>
  int jac_diffmat_grjm(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi Quadratures (point -1 included)
template<typename _Tp>
  _Tp jac_lagrange_grjm(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta);

/// Calculates the interpolation matrix for GRJM quadrature
template<typename _Tp>
  int jac_interpmat_grjm(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta);


/// Calculates the Gauss-Radau-Jacobi quadrature points (point +1 included)
template<typename _Tp>
  int jac_zeros_grjp(_Tp* z, const int Q, _Tp alpha, _Tp beta);

/// Calculates the Gauss-Radau-Jacobi quadrature weights (point +1 included)
template<typename _Tp>
  int jac_weights_grjp(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Gauss-Radau-Jacobi Derivative matrix (point +1 included)
template<typename _Tp>
  int jac_diffmat_grjp(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws);

/// Calculates the Lagrange polynomials through Gauss-Radau-Jacobi Quadratures (point +1 included)
template<typename _Tp>
  _Tp jac_lagrange_grjp(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta);

/// Calculates the interpolation matrix for GRJP quadrature
template<typename _Tp>
  int jac_interpmat_grjp(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta);

#endif /* __JACOBI_H__ */
