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
template<typename _Tp>
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


// Computes the jacobi polyniomials on several points
template<typename _Tp>
  int 
  jacobi_value_array(int np, const _Tp* x, _Tp* result_array,
		     int n, _Tp a, _Tp b)
  {
    if (n == 0)
      {
	for(int i = 0; i < np; ++i)
	  result_array[i] = _Tp{1};
	return 0;
      }

    if (n == 1)
      {
	for (int i = 0; i < np; ++i)
	  result_array[i] = 0.5 * (a - b + (a + b + _Tp{2}) * x[i]);
	return 0;
      }

    // General case,
    std::vector<_Tp> pnm1(np), pnm2(np, _Tp{1});
    for (int i = 0; i < np; ++i)
      pnm1[i] = 0.5 * (a - b + (a + b + _Tp{2}) * x[i]);

    // Start iterating:
    for (int k = 1; k < n; ++k)
      {
	auto a1 = _Tp{2} * (k + _Tp{1}) * (k + a + b + _Tp{1}) * (_Tp{2} * k + a + b);
	auto a2 = (_Tp{2} * k + a + b + _Tp{1}) * (a * a - b * b) / a1;
	auto a3 = (_Tp{2} * k + a + b) * (_Tp{2} * k + a + b + _Tp{1}) * (_Tp{2} * k + a + b + _Tp{2}) / a1;
	auto a4 = _Tp{2} * (k + a) * (k + b) * (_Tp{2} * k + a + b + _Tp{2}) / a1;
	for (int i = 0; i < np; ++i)
	  {
	    result_array[i] = (a2 + a3 * x[i]) * pnm1[i] - a4 * pnm2[i];
	    pnm2[i] = pnm1[i];
	    pnm1[i] = result_array[i];
	  }
      }

    return 0;
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
template<typename _Tp>
  int
  jacobi_deriv_array(int np, const _Tp* x, _Tp* result_array,
		     int n, double a, double b)
  {
    int code = jacobi_value_array(np, x, result_array,
			          n - 1, a + _Tp{1}, b + _Tp{1});
    if (code)
      return code;
    
    for (int i = 0; i < np; ++i)
      result_array[i] *= 0.5 * (a + b + n + _Tp{1});

    return code;
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
 * @param x A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::zeros_gj()
  {
    return this->jacobi_zeros(this->x.data(), this->Q, this->alpha, this->beta);
  }


/**
 * Calculates the Quadrature weights for Gauss-Jacobi integration.
 *
 * @param x Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code defined in gsl_errno.h
 */ 
template<typename _Tp>
  int
  jac_quadrature<_Tp>::weights_gj()
  {
    _Tp coef = std::pow(_Tp{2}, this->alpha + this->beta + _Tp{1})
	     * std::tgamma(this->alpha + _Tp(this->Q + 1))
	     / __gnu_cxx::factorial<_Tp>(this->Q)
	     * std::tgamma(this->beta + _Tp(this->Q + 1))
	     / std::tgamma(this->alpha + this->beta + _Tp(this->Q + 1));

    int code = jacobi_deriv_array(this->Q, this->x.data(), this->w.data(),
				  this->Q, this->alpha, this->beta);
    if (code)
      return code;

    for (int i = 0; i < this->Q; ++i)
      {
	const auto ww = this->w[i];
	const auto xx = this->x[i];
	this->w[i] = _Tp{1} / (ww * ww) * coef / (_Tp{1} - xx * xx);
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
 * @param x Quadrature nodes.
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code define in gsl_errno.h
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::diffmat_gj()
  {
    std::vector<_Tp> pnd(this->Q);

    jacobi_deriv_array(this->Q, this->x.data(), pnd.data(),
		       this->Q, this->alpha, this->beta);
    for (int i = 0; i < this->Q; ++i)
      for (int j = 0; j < this->Q; ++j)
	if (i != j)
	  this->D[i * this->Q + j] = (pnd[i] / pnd[j])
				   / (this->x[i] - this->x[j]);
	else
	  this->D[i * this->Q + j] = (this->alpha - this->beta
				      + (this->alpha + this->beta + _Tp{2}) * this->x[i])
				   / (_Tp{1} - this->x[i] * this->x[i]) / _Tp{2};

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
 * @param x The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_quadrature<_Tp>::lagrange_gj(int i, _Tp zz)
  {
    if (!compare(zz, this->x[i], EPS<_Tp>))
      return _Tp{1};

    return this->jacobi_value(zz, this->Q, this->alpha, this->beta)
	 / (this->jacobi_deriv(x[i], this->Q, this->alpha, this->beta) * (zz - this->x[i]));
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
 * @param imat A vector Q*xp.size() long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param x Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::interpmat_gj()
  {
    for (int i = 0; i < this->xp.size(); ++i)
      for (int j = 0; j < this->Q; ++j)
	this->imat[i * this->Q + j] = this->lagrange_gj(j, this->xp[i]);
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
 * @param x A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::zeros_glj()
  {
    // The zeros
    this->x[0] = _Tp{-1};
    this->x[Q - 1] = _Tp{1};

    return this->jacobi_zeros(this->x.data() + 1,
			      this->Q - 2, this->alpha + _Tp{1}, this->beta + _Tp{1});
  }


/**
 * Calculates the Quadrature weights for Gauss-Lobatto-Jacobi integration.
 *
 * @param x Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */ 
template<typename _Tp>
  int
  jac_quadrature<_Tp>::weights_glj()
  {
    _Tp coef = pow(_Tp{2}, this->alpha + this->beta + _Tp{1}) / _Tp(this->Q - 1)
	     * std::tgamma(this->alpha + this->Q)
	     / __gnu_cxx::factorial<_Tp>(this->Q - 1)
	     * std::tgamma(this->beta + this->Q)
	     / std::tgamma(this->alpha + this->beta + _Tp(this->Q + 1));

    jacobi_value_array(this->Q, this->x.data(), this->w.data(),
		       this->Q - 1, this->alpha, this->beta);

    this->w[0] = (this->beta + _Tp{1}) * coef / (this->w[0] * this->w[0]);
    for (int i = 1; i < this->Q - 1; ++i)
      {
	//const auto xx = this->x[i + 1];
	const auto ww = this->w[i + 1];
        this->w[i] = coef / (ww * ww);
      }
    this->w[this->Q - 1] = (this->alpha + _Tp{1}) * coef / std::pow(_Tp{2}, w[this->Q - 1]);

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
 * @param x Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::diffmat_glj()
  {
    std::vector<_Tp> pnd(this->Q);

    pnd[0] = pow(_Tp{-1}, this->Q)
	   * _Tp{2} * std::tgamma(this->Q + this->beta)
	   / std::tgamma(this->Q - 1)
	   / std::tgamma(this->beta + _Tp{2});

    pnd[this->Q - 1] = -_Tp{2} * std::tgamma(_Tp(this->Q) + this->alpha)
		     / std::tgamma(_Tp(this->Q - 1))
		     / std::tgamma(this->alpha + _Tp{2});
    jacobi_deriv_array(this->Q - 2, this->x.data() + 1, pnd.data() + 1,
		       this->Q - 2, this->alpha + _Tp{1}, this->beta + _Tp{1});

    for (int i = 1; i < this->Q - 1; ++i)
      pnd[i] *= (_Tp{1} - this->x[i]) * (_Tp{1} + this->x[i]);

    for (int i = 0; i < this->Q; ++i)
      for (int j = 0; j < this->Q; ++j)
	if (i != j)
	  this->D[i * this->Q + j] = (pnd[i] / pnd[j]) / (this->x[i] - this->x[j]);
	else
	  this->D[i * this->Q + i] = 0.5 * (this->alpha - this->beta + (this->alpha + this->beta) * this->x[i])
		       / (_Tp{1} - this->x[i] * this->x[i]);

    this->D[0] = 0.5 * (this->alpha - (this->Q - 1) * (this->Q + this->alpha + this->beta)) / (this->beta + _Tp{2});
    this->D[this->Q * this->Q - 1] = -0.5 * (this->beta - (this->Q - 1) * (this->Q + this->alpha + this->beta)) / (this->alpha + _Tp{2});

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
 * @param x The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_quadrature<_Tp>::lagrange_glj(int i, _Tp zz)
  {
    _Tp zi = this->x[i];
    if (!compare(zz, zi, EPS<_Tp>))
      return _Tp{1};

    return (_Tp{1} - zz * zz) * this->jacobi_value(zz, this->Q - 2, this->alpha + _Tp{1}, this->beta + _Tp{1})
	 / ((-2 * zi * this->jacobi_value(zi, this->Q - 2, this->alpha + _Tp{1}, this->beta + _Tp{1}) +
	    (_Tp{1} - zi * zi) * this->jacobi_deriv(zi, this->Q - 2, this->alpha + _Tp{1}, this->beta + _Tp{1})) * (zz - zi));
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
 * @param imat A vector Q*xp.size() long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param x Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::interpmat_glj()
  {
    for (int i = 0; i < this->xp.size(); ++i)
      for (int j = 0; j < this->Q; ++j)
	this->imat[i * this->Q + j] = this->lagrange_glj(j, this->xp[i]);
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
 * @param x A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::zeros_grjm()
  {
    // The zeros
    this->x[0] = _Tp{-1};

    return this->jacobi_zeros(this->x.data() + 1,
			      this->Q - 1, this->alpha, this->beta + _Tp{1});
  }

/**
 * Calculates the Quadrature weights for Gauss-Radau-Jacobi quadrature including point -1
 *
 * @param x Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */ 
template<typename _Tp>
  int
  jac_quadrature<_Tp>::weights_grjm()
  {
    _Tp coef = std::pow(_Tp{2}, this->alpha + this->beta)
	     / (this->beta + this->Q)
	     * std::tgamma(this->alpha + this->Q)
	     / __gnu_cxx::factorial<_Tp>(this->Q - 1)
	     * std::tgamma(this->beta + this->Q)
	     / std::tgamma(this->alpha + this->beta + _Tp(this->Q + 1));

    jacobi_value_array(this->Q, this->x.data(), this->w.data(),
		       this->Q - 1, this->alpha, this->beta);

    //this->w[0] = coef * (this->beta + _Tp{1}) * (_Tp{1} - xx) / (this->w[0] * this->w[0]);
    for (int i = 0/*1*/; i < this->Q; ++i)
      {
	const auto ww = this->w[i];
	const auto xx = this->x[i];
        w[i] = coef / (ww * ww) * (_Tp{1} - xx);
      }
    this->w[0] *= (this->beta + _Tp{1});

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
 * @param x Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::diffmat_grjm()
  {
    std::vector<_Tp> pnd(this->Q);

    jacobi_deriv_array(this->Q, this->x.data(), pnd.data(),
		       this->Q - 1, this->alpha, this->beta + _Tp{1});
    for (int i = 1; i < Q; ++i)
      pnd[i] *= (_Tp{1} + x[i]);

    pnd[0] = std::pow(_Tp{-1}, this->Q - 1)
	   * std::tgamma(this->beta + _Tp(this->Q + 1))
	   / std::tgamma(this->Q)
	   / std::tgamma(this->beta + _Tp{2});

    for (int i = 0; i < this->Q; ++i)
      for (int j = 0; j < this->Q; ++j)
	if (i != j)
	  this->D[i * this->Q + j] = (pnd[i] / pnd[j]) / (this->x[i] - this->x[j]);
	else
	  this->D[i * this->Q + j] = (this->alpha - this->beta + _Tp{1}
				   + (this->alpha + this->beta + _Tp{1}) * this->x[i])
				   / (_Tp{1} - this->x[i] * this->x[i]) / _Tp{2};
    this->D[0] = -0.5 * _Tp(this->Q - 1)
	       / (this->beta + _Tp{2})
	       * (this->alpha + this->beta + _Tp(this->Q + 1));

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
 * @param x The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_quadrature<_Tp>::lagrange_grjm(int i, _Tp zz)
  {
    _Tp zi = this->x[i];

    if (!compare(zz, this->x[i], EPS<_Tp>))
      return _Tp{1};

    return (_Tp{1} + zz) * this->jacobi_value(zz, this->Q - 1, this->alpha, this->beta + _Tp{1})
	 / ((this->jacobi_value(zi, this->Q - 1, this->alpha, this->beta + _Tp{1})
	    + (_Tp{1} + zi) * this->jacobi_deriv(zi, this->Q - 1, this->alpha, this->beta + _Tp{1})) * (zz - zi));
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
 * @param imat A vector Q*xp.size() long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param x Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::interpmat_grjm()
  {
    for (int i = 0; i < this->xp.size(); ++i)
      for (int j = 0; j < this->Q; ++j)
	this->imat[i * this->Q + j] = this->lagrange_grjm(j, this->xp[i]);
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
 * @param x A _Tp pointer to an array used to store the values of the zeros
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::zeros_grjp()
  {
    // The zeros
    this->x[Q - 1] = _Tp{1};

    return this->jacobi_zeros(this->x.data(),
			      this->Q - 1, this->alpha + _Tp{1}, this->beta);
  }


/**
 * Calculates the Quadrature weights for Gauss-Radau-Jacobi quadrature including point +1
 *
 * @param x Quadrature nodes
 * @param w Array containing the quadrature weights
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::weights_grjp()
  {
    _Tp coef = pow(_Tp{2}, this->alpha + this->beta)
	     / (this->alpha + this->Q)
	     * std::tgamma(this->alpha + this->Q)
	     / __gnu_cxx::factorial<_Tp>(this->Q - 1)
	     * std::tgamma(this->beta + this->Q)
	     / std::tgamma(this->alpha + this->beta + _Tp(this->Q + 1));

    jacobi_value_array(this->Q, this->x.data(), this->w.data(),
		       this->Q - 1, this->alpha, this->beta);
    for (int i = 0; i < this->Q; ++i)
      {
	const auto xx = this->x[i];
	const auto ww = this->w[i];
	w[i] = coef * (_Tp{1} + xx) / (ww * ww);
      }
    this->w[this->Q - 1] *= (this->alpha + _Tp{1});

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
 * @param x Quadrature nodes
 * @param D Derivative matrix, a vector with Q*Q _Tps
 * @param Q Number of quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::diffmat_grjp()
  {
    std::vector<_Tp> pnd(this->Q);

    jacobi_deriv_array(this->Q - 1, this->x.data(), pnd.data(),
		       this->Q - 1, this->alpha + _Tp{1}, this->beta);
    for (int i = 0; i < this->Q - 1; ++i)
      pnd[i] *= (_Tp{1} - x[i]);

    pnd[this->Q - 1] = -std::tgamma(this->alpha + _Tp(this->Q + 1))
		     / __gnu_cxx::factorial<_Tp>(this->Q - 1)
		     / std::tgamma(this->alpha + _Tp{2});
    
    for (int i = 0; i < this->Q; ++i)
      for (int j = 0; j < this->Q; ++j)
	if (i != j)
	  this->D[i * this->Q + j] = (pnd[i] / pnd[j]) / (x[i] - x[j]);
	else
	  this->D[i * this->Q + j] = (this->alpha - this->beta - _Tp{1}
				      + (this->alpha + this->beta + _Tp{1}) * this->x[i])
				   / (_Tp{1} - this->x[i] * this->x[i]) / _Tp{2};
    this->D[this->Q * this->Q - 1] = 0.5 * (this->Q - 1)
				   * (this->Q + this->alpha + this->beta + _Tp{1})
				   / (this->alpha + _Tp{2});
    
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
 * @param x The quadrature nodes
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return The value of the ith Lagrange polynomial at zz
 */
template<typename _Tp>
  _Tp
  jac_quadrature<_Tp>::lagrange_grjp(int i, _Tp zz)
  {
    auto zi = this->x[i];

    if (!compare(zz, this->x[i], EPS<_Tp>))
      return _Tp{1};

    return (_Tp{1} - zz) * this->jacobi_value(zz, this->Q - 1, this->alpha + _Tp{1}, this->beta)
			/ ((-this->jacobi_value(zi, this->Q - 1, this->alpha + _Tp{1}, this->beta)
	   + (_Tp{1} - zi) * this->jacobi_deriv(zi, this->Q - 1, this->alpha + _Tp{1}, this->beta)) * (zz - zi));
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
 * @param imat A vector Q*xp.size() long that will store the interpolation matrix
 * @param zp The points where the function should be interpolated
 * @param x Quadrature points
 * @param Q Number of quadrature points
 * @param alpha @f$\alpha@f$ parameter of Jacobi polynomial
 * @param beta @f$\beta@f$ parameter of Jacobi polynomial
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::interpmat_grjp()
  {
    for (int i = 0; i < this->xp.size(); ++i)
      for (int j = 0; j < this->Q; ++j)
	this->imat[i * this->Q + j] = this->lagrange_grjp(j, this->xp[i]);
    return 0;
  }


