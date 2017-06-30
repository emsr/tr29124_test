/**
 * gauss_jacobi_interface.tcc
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
#include <numeric>
#include "jacobi.h"
#include "integration_error.h"

/*
//  Y = alpha A X + beta Y
template<typename _Tp>
  void
  fake_gemv(// Specifies row-major (C) or column-major (Fortran) data ordering.
	    // SpecifiesIn file included from test_gauss_jacobi.cpp:14:0:
jacobi.h:111:10: note: declared private here
     _Tp* x;
          ^
test_gauss_jacobi.cpp:52:40: error: \u2018_Tp jac_quadrature<_Tp>::integrate(_Tp*) [with _Tp = double]\u2019 is private within this context
   auto integr = quad.integrate(f.data());
                                        ^
In file included from jacobi.h:238:0,
                 from test_gauss_jacobi.cpp:14:
gauss_jacobi_interface.tcc:291:3: note: declared private here
   jac_quadrature<_Tp>::integrate(_Tp* f)
   ^~~~~~~~~~~~~~~~~~~
gauss_jacobi_interface.tcc: In instantiation of \u2018_Tp jac_quadrature<_Tp>::integrate(_Tp*) [with _Tp = double]\u2019:
test_gauss_jacobi.cpp:52:40:   required from here
 whether to transpose matrix A.
	    int __M,	   // 
	    int __N,	   // 
	    _Tp __alpha,   // Scaling factor for the product of matrix A and vector X.
	    const _Tp*__A, // Matrix A.
	    int __lda,     // The size of the first dimension of matrix A; if you are passing a matrix A[m][n], the value should be m.
	    const _Tp*__X, // Vector X.
	    int __incX,    // Stride within X. For example, if incX is 7, every 7th element is used.
	    _Tp __beta,    // Scaling factor for vector Y.
	    _Tp* __Y,       // Vector Y
	    int __incY)    // Stride within Y. For example, if incY is 7, every 7th element is used.
  {
    if (__beta != _Tp{0})
      for (std::size_t __iy = 0, __ir = 0; __ir < __lda; __iy += __incY, ++__ir)
	__Y[__iy] *= __beta;

    if (__alpha != _Tp{0})
      for (std::size_t __iy = 0, __ir = 0; __ir < __lda; __iy += __incY, ++__ir)
	{
	  for (std::size_t __ix = 0, __ic = 0; __ic < __N; __ix += __incX, ++__ic)
	    __Y[__iy] += __alpha * __A[__ir][__ic] * __X[__ix];
	}
  }
*/
//    cblas_dgemv(CblasRowMajor, CblasNoTrans, quad->Q, quad->Q, _Tp{1}, quad->D, quad->Q, f, 1, _Tp{0}, d, 1);
//    matvec(quad->Q, quad->D, f, d);

//  Y = A X
template<typename _Tp>
  void
  matvec(std::size_t __n, const _Tp* __A, const _Tp* __x, _Tp* __y)
  {
    for (std::size_t __ir = 0; __ir < __n; ++__ir)
      {
	__y[__ir] = _Tp{0};
	for (std::size_t __ic = 0; __ic < __n; ++__ic)
	  __y[__ir] += __A[__ir][__ic] * __x[__ic];
      }
  }


/**
 * This function calculates the quadrature nodes weights and derivative matrix
 * and stores the data on a previously allocated jac_quadrature structure.
 *
 * @param quad A previously allocated jac_quadrature structure
 * @param qtype Quadrature type
 * @param a Alpha weight
 * @param b Beta weight
 * @return An error code or 0
 */
template<typename _Tp>
  int
  jac_quadrature<_Tp>::quadrature_zwd()
  {
    // Calculates the zeros of the quadrature
    int err = 0;
    switch (this->type)
      {
      case Gauss:
	err = this->zeros_gj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = this->weights_gj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = this->diffmat_gj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Lobatto:
	err = this->zeros_glj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = this->weights_glj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = this->diffmat_glj();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Radau_lower:
	err = this->zeros_grjm();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = this->weights_grjm();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = this->diffmat_grjm();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      case Gauss_Radau_upper:
	err = this->zeros_grjp();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the zeros", err, _Tp{}, _Tp{});
	err = this->weights_grjp();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the weightd", err, _Tp{}, _Tp{});
	err = this->diffmat_grjp();
	if (err)
	  __gnu_test::__throw__IntegrationError("Problem calculating the differentiation matrix", err, _Tp{}, _Tp{});
	break;
      default:
	__gnu_test::__throw__IntegrationError("Illegal quadrature type", err, _Tp{}, _Tp{});
      }

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
  jac_quadrature<_Tp>::interpmat_alloc(int npoints, _Tp* xp)
  {
    if (npoints < 1)
      std::__throw_domain_error("The number of interpolating points should be at least 1");

    this->xp.resize(npoints);
    this->imat.resize(this->Q * npoints);
    for (int i = 0; i < npoints; ++i)
      this->xp[i] = xp[i];

    int err = 0;

    switch (this->type)
    {
    case Gauss:
      err = this->interpmat_gj();
      break;
    case Gauss_Lobatto:
      err = this->interpmat_glj();
      break;
    case Gauss_Radau_lower:
      err = this->interpmat_grjm();
      break;
    case Gauss_Radau_upper:
      err = this->interpmat_grjp();
      break;
    default:
      __gnu_test::__throw__IntegrationError("Illegal quadrature type", err, _Tp{}, _Tp{});
    }

    return err;
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
  template<typename _Func>
    _Tp
    jac_quadrature<_Tp>::integrate(_Func fun)
    {
      std::vector<_Tp> f(this->Q);
      for (int i = 0; i < this->Q; ++i)
        f[i] = fun(this->x[i]);
      return std::inner_product(std::begin(this->w), std::end(this->w),
				std::begin(f), _Tp{0});
    }
    

/**
 * Calculates the derivative of a function known at quadrature points given the derivative matrix
 * The derivative matrix should have been calculated before using one of the functions *_diffmat.
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
  jac_quadrature<_Tp>::differentiate(_Tp* f, _Tp* d)
  {
    matvec(this->Q, this->D, f, d);
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
  jac_quadrature<_Tp>::interpolate(_Tp* f, _Tp* fout)
  {
    if (this->xp.size() == 0)
      std::__throw_runtime_error("No interpolation info was setup");

    //cblas_dgemv(CblasRowMajor, CblasNoTrans, this->np, this->Q, _Tp{1}, this->imat, this->Q, f, 1, _Tp{0}, fout, 1);
    matvec(this->Q, this->D, f, fout);
    return 0;
  }
