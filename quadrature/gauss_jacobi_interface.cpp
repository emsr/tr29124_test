/**
 * jacobi.c
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

/** @file  Wrappers to basic functions defined in gauss_quad.c
 */

#include <cstdlib>
#include "jacobi.h"
#include "integration_error.h"


/**
 * Allocates memory for a jac_quadrature structure.
 * @param qtype Quadrature type
 * @param nq Number of quadrature nodes
 * @param a Alpha weight
 * @param b Beta weight
 * @return A pointer to newly created jac_quadrature structure. nullptr if some problem occurred
 */
template<typename _Tp>
  jac_quadrature<_Tp>*
  jac_quadrature_alloc(int nq)
  {
    if (nq < 1)
      std::__throw_domain_error("The number of quadrature points should be at least 1");

    void* mem_block = malloc(sizeof(jac_quadrature<_Tp>) + (2 * nq + nq * nq) * sizeof(_Tp));

    if (mem_block == nullptr)
      std::__throw_runtime_error("Memory allocation was not possible");

    jac_quadrature<_Tp>* quad = static_cast<jac_quadrature<_Tp>*>(mem_block);

    quad->Q = nq;

    quad->x = static_cast<_Tp*>(mem_block + sizeof(jac_quadrature<_Tp>));
    quad->w = quad->x + nq;
    quad->D = quad->x + 2 * nq;
    quad->np = 0;
    quad->Imat = nullptr;
    quad->xp = nullptr;

    return quad;
  }

/**
 * Releases memory allocated by jac_quadrature_alloc. It also releases memory allocated by the
 * jac_interpmat_alloc if necessary
 *
 * @param quad A jac_quadrature structure allocated with jac_quadrature_alloc
 */
template<typename _Tp>
  void 
  jac_quadrature_free(jac_quadrature<_Tp>* quad)
  {
    jac_interpmat_free(quad);
    free(quad);
  }


/**
 * This function calculates the quadrature nodes weights and derivative matrix
 * and stores the data on a previously allocated jac_quadrature structure.
 *
 * @param quad A previously allocated jac_quadrature structure
 * @param qtype Quadrature type
 * @param a Alpha weight
 * @param b Beta weight
 * @param ws Workspace 3*quad->Q long used in calculations. If it is nullptr, memory is allocated and at the end released
 * @return An error code or 0
 * 
 */
template<typename _Tp>
  int 
  jac_quadrature_zwd(jac_quadrature<_Tp>* quad, enum jac_quad_type qtype, _Tp a, _Tp b, _Tp* ws)
  {
    int allocated = 0;
    if (ws == nullptr)
      {
	// Try to allocate the workspace:
	ws = static_cast<_Tp*>(malloc(3 * quad->Q));
	if (ws == nullptr)
	  std::__throw_runtime_error("Could not allocate workspace memory");
	allocated = 1;
      }
    quad->alpha = a;
    quad->beta = b;
    quad->type = qtype;
    int err = 0;

    // Calculates the zeros of the quadrature
    switch (quad->type)
      {
      case Gauss:
	err = jac_zeros_gj(quad->x, quad->Q, quad->alpha, quad->beta);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = jac_weights_gj(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = jac_diffmat_gj(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Lobatto:
	err = jac_zeros_glj(quad->x, quad->Q, quad->alpha, quad->beta);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = jac_weights_glj(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = jac_diffmat_glj(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Radau_m1:
	err = jac_zeros_grjm(quad->x, quad->Q, quad->alpha, quad->beta);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = jac_weights_grjm(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = jac_diffmat_grjm(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Radau_p1:
	err = jac_zeros_grjp(quad->x, quad->Q, quad->alpha, quad->beta);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = jac_weights_grjp(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = jac_diffmat_grjp(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      default:
	__gnu_test::__throw__IntegrationError("Illegal quadrature type", err, _Tp{}, _Tp{});
      }

    if (allocated)
      free(ws);

    return 0;
  }


/**
 * Allocates memory for interpolation matrix for a jac_quadrature structure. 
 * It also copies the interpolation points.
 *
 * @param quad A previously allocated jac_quadrature structure
 * @param npoints Number of interpolation points
 * @param xp Interpolation points
 * @return Error code or 0
 */
template<typename _Tp>
  int 
  jac_interpmat_alloc(jac_quadrature<_Tp>* quad, int npoints, _Tp* xp)
  {
    if (npoints < 1)
      std::__throw_domain_error("The number of interpolating points should be at least 1");

    // Allocate memory for the interpolation matrix and points
    quad->xp = static_cast<_Tp*>(malloc(sizeof(_Tp) * (npoints + npoints*quad->Q)));
    if (!quad->xp)
      std::__throw_runtime_error("Memory for interpolation matrix could not be allocated");

    quad->Imat = quad->xp + npoints;
    for (int i = 0; i < npoints; ++i)
      quad->xp[i] = xp[i];

    int err = 0;

    quad->np = npoints;

    switch (quad->type)
    {
    case Gauss:
      err = jac_interpmat_gj(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case Gauss_Lobatto:
      err = jac_interpmat_glj(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case Gauss_Radau_m1:
      err = jac_interpmat_grjm(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case Gauss_Radau_p1:
      err = jac_interpmat_grjp(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    default:
      __gnu_test::__throw__IntegrationError("Illegal quadrature type", err, _Tp{}, _Tp{});
    }

    return err;
  }


/**
 * Frees memory used by jac_ionterpmat_alloc
 *
 * @param quad An allocated jac_quadrature structure.
 */
template<typename _Tp>
  void
  jac_interpmat_free(jac_quadrature<_Tp>* quad)
  {
    if (quad->xp)
      free(quad->xp);
    quad->np = 0;
    quad->xp = nullptr;
    quad->Imat = nullptr;
  }


/** Aproximates the integral of a function using a a quadrature.
 *
 * @f[
 *  \int_{-1}^1 f(x) dx \approx \sum_{i=0}^{Q-1} w_i f(x_i)
 * @f]
 *
 * Where @f$x_i@f$ and @f$w_i@f$ are given by the quadrature method emplyed.
 *
 * @param quad A structure containing quadrature information
 * @param f Value of function at quadrature points
 * @return @f$\int_{-1}^{1} f(x) dx \approx \sum_{i=0}^{Q-1} w_i f(x_i)@f$
 */
template<typename _Tp>
  _Tp 
  jac_integrate(jac_quadrature<_Tp>* quad, _Tp* f)
  {
    auto sum = cblas_ddot(quad->Q, quad->w, 1, f, 1);
    return sum;
  }
    

/**
 * Calculates the derivative of a function known at quadrature points given the derivative matrix
 * The derivative matrix should have been calculated before using one of the functions jac_**_diffmat.
 *
 * The derivative is calculated according to the following equation:
 * 
 * @f[
 *    \left.\frac{du(x)}{dx}\right|_{x=x_i} = \sum_{j=0}^{Q-1} D_{ij} u(x_j)
 * @f]
 *
 * In this equation,
 * @f[
 *   D_{ij} = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
 * @f]
 *
 *  where @f$h_j(x)@f$ is the Lagrange polynomial through the jth quadrature node
 *  
 * @param quad A structure containing quadrature information
 * @param f Value of function at quadrature points
 * @param d The Estimated derivative at quadrature points: 
 * @return 0 if everything was ok. Otherwise return an error code
 */
template<typename _Tp>
  int
  jac_differentiate(jac_quadrature<_Tp>* quad, _Tp* f, _Tp* d)
  {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, quad->Q, quad->Q, _Tp{1}, quad->D, quad->Q, f, 1, _Tp{0}, d, 1);
    return 0;
  }


/**
 * Interpolates a function given at quadrature points using the interpolation matrix.
 * The interpolation matrix is given by the following equation:
 * @f[
 *   I_{ij} = h_j(x_i)
 * @f]
 *
 * where @f$x_i@f$ are the points where the function should be interpolated.
 * With this matrix, the interpolation consists of the following operation:
 * @f[
 *    f^{interpolated}(x_i) = \sum_{j=0}^{Q-1} f_j h_j(x_i)
 * @f]
 *
 * @param quad Quadrature information (including interpolation matrix)
 * @param f Value of function at quadrature points
 * @param fout Interpolated values
 * @return 0 if everything was ok. Otherwise return an error code
 * 
 */
template<typename _Tp>
  int
  jac_interpolate(jac_quadrature<_Tp>* quad, _Tp* f, _Tp* fout)
  {
    if (!quad->np)
      std::__throw_runtime_error("No interpolation info was setup");

    cblas_dgemv(CblasRowMajor, CblasNoTrans, quad->np, quad->Q, _Tp{1}, quad->Imat, quad->Q, f, 1, _Tp{0}, fout, 1);
    return 0;
  }
