/* gauss_quad.c
 * 
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

/**
  @file gauss_jacobi_integrate.cpp
  @brief This file implements Gauss-Jacobi quadrature related functions

  This library implements 4 types of quadrature
*/
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <limits>

#include "jacobi.h"
#include "integration_error.h"

/// Admissible convergence error
const _Tp EPS = _Tp(300) * std::numeric_limits<_Tp>::epsilon();


template<typename _Tp>
  static int
  compare(_Tp x, _Tp y, _Tp eps)
  {
    const auto dx = x - y;

    if (std::abs(dx) < eps)
      return 0;

    return dx > 0 ? 1 : -1;
  }


// **** GAUSS QUADRATURE ****


/**
 * Finds the zeros of Gauss-Jacobi quadrature
 *
 * It computes the following equation
 * @f[
 *   P^{(\alpha, \beta)}_Q(x_i) = 0
 * @f]
 *
 * @param z A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_zeros_gj(_Tp* z, const int Q, _Tp alpha, _Tp beta)
  {
    return jac_jacobi_zeros(z, Q, alpha, beta);
  }


/**
 * Calculates the Quadrature weights for Gauss-Jacobi integration.
 *
 * @param z Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with room for 2*Q _Tps. If it is null, the memory will be allocated with malloc and at the end released
 * @return 0 if everything was ok. Otherwise return an error code defined in gsl_errno.h
 */ 
template<typename _Tp>
  int
  jac_weights_gj(_Tp* z, _Tp* w, const int Q, const _Tp alpha, const _Tp beta, _Tp* ws)
  {
    _Tp coef = std::pow(_Tp{2}, alpha + beta + _Tp{1})
	     * (std::tgamma(alpha + _Tp{Q + 1}) / factorial(Q))
	     * (std::tgamma(beta + _Tp{Q + 1}) / std::tgamma(alpha + beta + _Tp{Q + 1}));

    int code = jac_djacobi_array(Q, z, Q, w, alpha, beta, ws);
    if (code)
      return code;

    for (int i = 0; i < Q; ++i)
      {
	const auto ww = w[i];
	const auto x = z[i];
	w[i] = _Tp{1} / (ww * ww) * coef / (_Tp{1} - x * x);
      }

    return 0;
  }


/**
 * Calculates the derivative matrix for Gauss-Jacobi quadrature.
 * The matrix should be preallocated, using a vector with Q*Q _Tps
 * This function calculates
 * @f[
 *   D[j+iQ] = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
 * @f]
 *
 * @param z Quadrature nodes.
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with at least 3*Q _Tps. If it is null, the memory will be allocated with malloc and at the end released.
 * @return 0 if everything was ok. Otherwise return an error code define in gsl_errno.h
 */
template<typename _Tp>
  int
  jac_diffmat_gj(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp* jac_deriv;
    _Tp* pnm1, *pnm2;
    bool mem_allocated = false;
    if (ws == nullptr)
      {
	jac_deriv = static_cast<_Tp*>(malloc(3 * Q * sizeof(_Tp)));
	if (!jac_deriv)
	  return GSL_ENOMEM;

	mem_allocated = true;
	pnm1 = jac_deriv + Q;
	pnm2 = pnm1 + Q;
      }
    else
      {
	jac_deriv = ws;
	pnm1 = ws+Q;
	pnm2 = pnm1+Q;
      }


    jac_djacobi_array(Q, z, Q, jac_deriv, alpha, beta, pnm1);
    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < Q; ++j)
	if (i != j)
	  D[i * Q + j] = (jac_deriv[i] / jac_deriv[j]) / (z[i] - z[j]);
	else
	  D[i * Q + j] = (alpha - beta + (alpha + beta + _Tp{2}) * z[i])
		       / (_Tp{1} - z[i] * z[i]) / _Tp{2};

    if (mem_allocated)
      free(jac_deriv);

    return 0;
  }


/**
 * This function calculates the ith Gauss-Jacobi lagrange interpolants at a point zz:
 * @f[
 *   h_i(zz)
 * @f]
 *
 * @param i The reference node of the Lagrange polynomial
 * @param zz The point where the Lagrange polynomial should be calculated
 * @param Q Number of quadrature nodes
 * @param z The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_lagrange_gj(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta)
  {
    if (!compare(zz, z[i], EPS))
      return _Tp{1};

    return jac_jacobi(zz, Q, alpha, beta)
	 / (jac_djacobi(z[i], Q, alpha, beta) * (zz - z[i]));
  }


/**
 * This function calculates the interpolation matrix for Gauss-Jacobi quadrature.
 * The interpolation matrix is simply the value of the Lagrange polynomials for each
 * point where the function should be interpolated, for each Lagrange polynomial:
 *
 * @f[
 *    I_{ij} = I[iQ+j] = h_j(x_i)
 * @f]
 *
 * @param imat A vector np*Q long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param np Number of points where the function should be interpolated
 * @param z Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_interpmat_gj(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta)
  {
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < Q; ++j)
	imat[i * Q + j] = jac_lagrange_gj(j, zp[i], Q, z, alpha, beta);
    return 0;
  }


// **** GAUSS-LOBATTO QUADRATURE ****


/**
 * Finds the zeros of Gauss-Lobatto-Jacobi quadrature
 *
 * This quadrature includes the endpoints (-1, 1). The other quadrature nodes are given by the zeros of
 * @f[
 * P^{\alpha+1, \beta+1}_{Q-2}(x_i) = 0
 * @f]
 *
 * @param z A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_zeros_glj(_Tp* z, const int Q, _Tp alpha, _Tp beta)
  {
    // The zeros
    z[0] = -_Tp{1};
    z[Q - 1] = _Tp{1};

    return jac_jacobi_zeros(z + 1, Q - 2, alpha + _Tp{1}, beta + _Tp{1});
  }


/**
 * Calculates the Quadrature weights for Gauss-Lobatto-Jacobi integration.
 *
 * @param z Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 2*Q _Tps. If it is null, the memory will be allocated with malloc
 * @return 0 if everything was ok. Otherwise return an error code
 */ 
template<typename _Tp>
  int
  jac_weights_glj(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp coef = pow(_Tp{2}, alpha + beta + _Tp{1}) / (Q - 1)
	     * (std::tgamma(alpha+Q) / gsl_sf_fact(Q - 1) )
	     * (std::tgamma(beta + Q) / std::tgamma(alpha + beta + Q + _Tp{1}));
    _Tp* ww, *x;

    jac_jacobi_array(Q, z, Q - 1, w, alpha, beta, ws);
    w[0] = (beta + _Tp{1}) * coef / (w[0] * w[0]);
    w[Q - 1] = (alpha + _Tp{1}) * coef / std::pow(_Tp{2}, w[Q - 1]);
    for (int i = 1, ww = w + 1, x = z + 1; i < Q - 1; ++i, ++ww, ++x)
      *ww = coef / (*ww * *ww);

    return 0;
  }


/**
 * Calculates the derivative matrix for Gauss-Lobatto-Jacobi quadrature.
 * The matrix should be preallocated, using a vector with Q*Q _Tps
 * This function calculates
 * @f[
 *   D[j+iQ] = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
 * @f]
 *
 * @param z Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 3*Q _Tps. If it is null, the memory will be allocated with malloc
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_diffmat_glj(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp* pnm1, *pnm2, *pqd;
    bool mem_allocated = false;
    if (ws == nullptr)
      {
	pqd = static_cast<_Tp*>(malloc(3 * Q * sizeof(_Tp)));
	if (!pqd)
	  return GSL_ENOMEM;

	mem_allocated = true
	pnm1 = pqd + Q;
	pnm2 = pnm1 + Q;
      }
    else
      {
	pqd = ws;
	pnm1 = ws + Q;
	pnm2 = pnm1 + Q;
      }

    pqd[0] = pow(-_Tp{1}, Q) * _Tp{2} * std::tgamma(Q + beta)
	   / (std::tgamma(Q-1) * std::tgamma(beta + _Tp{2}));

    pqd[Q-1] = -_Tp{2} * std::tgamma(Q+alpha)
	     / (std::tgamma(Q-1)*std::tgamma(alpha + _Tp{2}));
    jac_djacobi_array(Q - 2, z + 1, Q - 2, pqd + 1, alpha + _Tp{1}, beta + _Tp{1}, pnm1);

    for (int i = 1; i < Q - 1; ++i)
      pqd[i] *= (_Tp{1} - z[i]) * (_Tp{1} + z[i]);

    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < Q; ++j)
	if (i != j)
	  D[i * Q + j] = (pqd[i] / pqd[j]) / (z[i] - z[j]);
	else
	  D[i * Q + i] = 0.5 * (alpha - beta + (alpha + beta) * z[i])
		       / (_Tp{1} - z[i] * z[i]);

    D[0] = 0.5 * (alpha - (Q - 1) * (Q + alpha + beta)) / (beta + _Tp{2});
    D[Q * Q - 1] = -0.5 * (beta - (Q - 1) * (Q + alpha + beta)) / (alpha + 2);

    if (mem_allocated)
      free(pqd);

    return 0;
  }


/**
 * This function calculates the ith Gauss-Lobatto-Jacobi lagrange interpolants at a point zz:
 * @f[
 *   h_i(zz)
 * @f]
 *
 * @param i The reference node of the Lagrange polynomial
 * @param zz The point where the Lagrange polynomial should be calculated
 * @param Q Number of quadrature nodes
 * @param z The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp 
  jac_lagrange_glj(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta)
  {
    _Tp zi = z[i];
    if (!compare(zz, zi, EPS))
      return _Tp{1};

    return (_Tp{1} - zz * zz) * jac_jacobi(zz, Q - 2, alpha + _Tp{1}, beta + _Tp{1})
	 / ((-2 * zi * jac_jacobi(zi, Q - 2, alpha + _Tp{1}, beta + _Tp{1}) +
	    (_Tp{1} - zi * zi) * jac_djacobi(zi, Q - 2, alpha + _Tp{1}, beta + _Tp{1})) * (zz - zi));
  }


/**
 * This function calculates the interpolation matrix for Gauss-Lobatto-Jacobi quadrature.
 * The interpolation matrix is simply the value of the Lagrange polynomials for each
 * point where the function should be interpolated, for each Lagrange polynomial:
 *
 * @f[
 *    I_{ij} = I[iQ+j] = h_j(x_i)
 * @f]
 *
 * @param imat A vector np*Q long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param np Number of points where the function should be interpolated
 * @param z Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_interpmat_glj(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta)
  {
    for (int i = 0; i < np; ++i)
      for (jint  = 0; j < Q; ++j)
	imat[i * Q + j] = jac_lagrange_glj(j, zp[i], Q, z, alpha, beta);
    return 0;
  }


// **** GAUSS-RADAU -1 QUADRATURE ****


/**
 * Finds the zeros of Gauss-Radau-Jacobi quadrature
 *
 * This quadrature includes the endpoint -1. The other quadrature nodes are given by the zeros of
 * @f[
 * P^{(\alpha, \beta+1)}_{Q-1}(x_i) = 0
 * @f]
 *
 * @param z A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_zeros_grjm(_Tp* z, const int Q, _Tp alpha, _Tp beta)
  {
    // The zeros
    z[0] = -_Tp{1};

    return jac_jacobi_zeros(z + 1, Q - 1, alpha, beta + _Tp{1});
  }

/**
 * Calculates the Quadrature weights for Gauss-Radau-Jacobi quadrature including point -1
 *
 * @param z Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 2*Q _Tps. If it is null, the memory will be allocated with malloc and at the ens released
 * @return 0 if everything was ok. Otherwise return an error code
 */ 
template<typename _Tp>
  int 
  jac_weights_grjm(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp coef = pow(_Tp{2}, alpha + beta)
	     / (beta + Q) * (std::tgamma(alpha + Q) / factorial(Q - 1) )
	     * (std::tgamma(beta + Q) / std::tgamma(alpha + beta + Q + _Tp{1}));
    _Tp* ww, *x;

    jac_jacobi_array(Q, z, Q - 1, w, alpha, beta, ws);
    for (int i = 0, ww = w, x = z; i < Q; ++i, ++ww, ++x)
      *ww = coef / (*ww * *ww) * (1 - *x);

    w[0] *= (beta + _Tp{1});

    return 0;
  }


/**
 * Calculates the derivative matrix for Gauss-Radau-Jacobi quadrature including the end point -1.
 * The matrix should be preallocated, using a vector with Q*Q _Tps
 * This function calculates
 * @f[
 *   D[j+iQ] = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
 * @f]
 *
 * @param z Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 3*Q _Tps. If it is null, the memory will be allocated with malloc
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_diffmat_grjm(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp* pnm1, *pnm2, *pqd;
    bool mem_allocated = false;
    if (ws == nullptr)
      {
	pqd = static_cast<_Tp*>(malloc(3 * Q * sizeof(_Tp)));
	if (!pqd)
	  std::__throw_runtime_error("no memory");

	mem_allocated = true;
	pnm1 = pqd + Q;
	pnm2 = pnm1 + Q;
      }
    else
      {
	pqd = ws;
	pnm1 = pqd + Q;
	pnm2 = pnm1 + Q;
      }

    jac_djacobi_array(Q, z, Q - 1, pqd, alpha, beta + _Tp{1}, pnm1);
    for (int i = 1; i < Q; ++i)
      pqd[i] *= (_Tp{1} + z[i]);

    pqd[0] = std::pow(-_Tp{1}, Q - 1) * std::tgamma(beta + _Tp{Q + 1})
	   / std::tgamma(Q) / std::tgamma(beta + _Tp{2});

    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < Q; ++j)
	if (i != j)
	  D[i * Q + j] = (pqd[i] / pqd[j]) / (z[i] - z[j]);
	else
	  D[i * Q + j] = (alpha - beta + _Tp{1} + (alpha + beta + _Tp{1}) * z[i])
		       / (_Tp{1} - z[i] * z[i]) / _Tp{2};
    D[0] = -0.5 * (Q - 1) / (beta + _Tp{2}) * (alpha + beta + _Tp{Q + 1});

    if (mem_allocated)
      free(pqd);

    return 0;
  }


/**
 * This function calculates the ith Gauss-Radau-Jacobi (including end point -1)
 * Lagrange interpolants at a point zz:
 * @f[
 *   h_i(zz)
 * @f]
 *
 * @param i The reference node of the Lagrange polynomial
 * @param zz The point where the Lagrange polynomial should be calculated
 * @param Q Number of quadrature nodes
 * @param z The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp 
  jac_lagrange_grjm(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta)
  {
    _Tp zi = z[i];

    if (!compare(zz, z[i], EPS)) return _Tp{1};

    return (_Tp{1} + zz) * jac_jacobi(zz, Q - 1, alpha, beta + _Tp{1})
	 / ((jac_jacobi(zi, Q - 1, alpha, beta + _Tp{1})
	    + (_Tp{1} + zi) * jac_djacobi(zi, Q - 1, alpha, beta + _Tp{1})) * (zz - zi));
  }


/**
 * This function calculates the interpolation matrix for Gauss-Radau-Jacobi (including end point -1) quadrature.
 * The interpolation matrix is simply the value of the Lagrange polynomials for each
 * point where the function should be interpolated, for each Lagrange polynomial:
 *
 * @f[
 *    I_{ij} = I[iQ+j] = h_j(x_i)
 * @f]
 *
 * @param imat A vector np*Q long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param np Number of points where the function should be interpolated
 * @param z Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_interpmat_grjm(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta)
  {
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < Q; ++j)
	imat[i * Q + j] = jac_lagrange_grjm(j, zp[i], Q, z, alpha, beta);
    return 0;
  }


// **** GAUSS-RADAU +1 QUADRATURE ****


/**
 * Finds the zeros of Gauss-Radau-Jacobi quadrature
 *
 * This quadrature includes the endpoint 1.
 * The other quadrature nodes are given by the zeros of
 * @f[
 *   P^{(\alpha+1, \beta)}_{Q-1}(x_i) = 0
 * @f]
 *
 * @param z A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int 
  jac_zeros_grjp(_Tp* z, const int Q, _Tp alpha, _Tp beta)
  {
    // The zeros
    z[Q - 1] = _Tp{1};

    return jac_jacobi_zeros(z, Q - 1, alpha + _Tp{1}, beta);
  }


/**
 * Calculates the Quadrature weights for Gauss-Radau-Jacobi quadrature including point +1
 *
 * @param z Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 2*Q _Tps. If it is null, the memory will be allocated with malloc
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_weights_grjp(_Tp* z, _Tp* w, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp coef = pow(_Tp{2}, alpha + beta)
	     / (alpha + Q) * (std::tgamma(alpha + Q) / factorial(Q - 1))
	     * (std::tgamma(beta + Q) / std::tgamma(alpha + beta + Q + _Tp{1}));
    _Tp* ww, *x;

    jac_jacobi_array(Q, z, Q - 1, w, alpha, beta, ws);
    for (int i = 0, ww = w, x = z; i < Q; ++i, ++ww, ++x)
      *ww = coef / (*ww * *ww) * (1 + *x);

    w[Q - 1] *= (alpha + _Tp{1});

    return 0;
  }


/**
 * Calculates the derivative matrix for Gauss-Radau-Jacobi quadrature including the end point +1.
 * The matrix should be preallocated, using a vector with Q * Q _Tps
 * This function calculates
 * @f[
 *   D[j+iQ] = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
 * @f]
 *
 * @param z Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @param ws Workspace with 3 * Q _Tps. If it is null, the memory will be allocated with malloc
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_diffmat_grjp(_Tp* z, _Tp* D, const int Q, _Tp alpha, _Tp beta, _Tp* ws)
  {
    _Tp* pnm1, *pnm2, *pqd;
    bool mem_allocated = false;
    if (ws == nullptr)
      {
	pqd = static_cast<_Tp*>()malloc(3 * Q * sizeof(_Tp));
	if (!pqd)
	  return GSL_ENOMEM;
	
	mem_allocated = true;
	pnm1 = pqd + Q;
	pnm2 = pnm1 + Q;
      }
    else
      {
	pqd = ws;
	pnm1 = pqd + Q;
	pnm2 = pnm1 + Q;
      }

    jac_djacobi_array(Q - 1, z, Q-1, pqd, alpha + _Tp{1}, beta, pnm1);
    for (int i = 0; i < Q - 1; ++i)
	pqd[i] *= (_Tp{1} - z[i]);

    pqd[Q-1] = -std::tgamma(alpha + _Tp{Q + 1})
	     / factorial(Q - 1) / std::tgamma(alpha + _Tp{2});
    
    for (int i = 0; i < Q; ++i)
      for (int j = 0; j < Q; ++j)
	if (i != j)
	  D[i * Q + j] = (pqd[i] / pqd[j]) / (z[i] - z[j]);
	else
	  D[i * Q + j] = (alpha - beta - _Tp{1} + (alpha + beta + _Tp{1}) * z[i])
		       / (_Tp{1} - z[i] * z[i]) / _Tp{2};
    D[Q * Q - 1] = 0.5 * (Q - 1) * (Q + alpha + beta + _Tp{1}) / (alpha + _Tp{2});

    if (mem_allocated)
      free(pqd);
    
    return 0;
  }


/**
 * This function calculates the ith Gauss-Radau-Jacobi (including end point +1)
 * Lagrange interpolants at a point zz:
 * @f[
 *   h_i(zz)
 * @f]
 *
 * @param i The reference node of the Lagrange polynomial
 * @param zz The point where the Lagrange polynomial should be calculated
 * @param Q Number of quadrature nodes
 * @param z The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_lagrange_grjp(int i, _Tp zz, int Q, _Tp* z, _Tp alpha, _Tp beta)
  {
    auto zi = z[i];

    if (!compare(zz, z[i], EPS))
      return _Tp{1};

    return (_Tp{1} - zz) * jac_jacobi(zz, Q - 1, alpha + _Tp{1}, beta)
	 / ((-jac_jacobi(zi, Q - 1, alpha + _Tp{1}, beta)
	   + (_Tp{1} - zi) * jac_djacobi(zi, Q - 1, alpha + _Tp{1}, beta)) * (zz - zi));
  }


/**
 * This function calculates the interpolation matrix for Gauss-Radau-Jacobi (including end point +1) quadrature.
 * The interpolation matrix is simply the value of the Lagrange polynomials for each
 * point where the function should be interpolated, for each Lagrange polynomial:
 *
 * @f[
 *    I_{ij} = I[iQ+j] = h_j(x_i)
 * @f]
 *
 * @param imat A vector np*Q long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param np Number of points where the function should be interpolated
 * @param z Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_interpmat_grjp(_Tp* imat, _Tp* zp, int np, _Tp* z, int Q, _Tp alpha, _Tp beta)
  {
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < Q; ++j)
	imat[i * Q + j] = jac_lagrange_grjp(j, zp[i], Q, z, alpha, beta);
    return 0;
  }


