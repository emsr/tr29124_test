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

/** \file  Wrappers to basic functions defined in gauss_quad.c
 */

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include "jacobi.h"


/**
   Allocates memory for a jac_quadrature structure.
   \param qtype Quadrature type
   \param nq Number of quadrature nodes
   \param a Alpha weight
   \param b Beta weight
   \return A pointer to newly created jac_quadrature structure. NULL if some problem occurred
 */
jac_quadrature *
jac_quadrature_alloc(int nq)
{
  
  if (nq < 1)
    {
      GSL_ERROR_NULL("The number of quadrature points should be at least 1", GSL_EINVAL);
    }


  void *mem_block = malloc(sizeof(jac_quadrature) + (2*nq+nq*nq)*sizeof(double));

  if (mem_block == NULL)
    {
      GSL_ERROR_NULL("Memory allocation was not possible", GSL_ENOMEM);
    }

  jac_quadrature *quad = (jac_quadrature *) mem_block;

  quad->Q = nq;

  quad->x = (double *) (mem_block + sizeof(jac_quadrature));
  quad->w = quad->x + nq;
  quad->D = quad->x + 2*nq;
  quad->np = 0;
  quad->Imat = NULL;
  quad->xp = NULL;
  
  
  return quad;
  
  
}

/**
  Releases memory allocated by jac_quadrature_alloc. It also releases memory allocated by the
  jac_interpmat_alloc if necessary

  \param quad A jac_quadrature structure allocated with jac_quadrature_alloc
 */
void 
jac_quadrature_free(jac_quadrature *quad)
{

  jac_interpmat_free(quad);
  free(quad);
  
}


  
  
/**
   This function calculates the quadrature nodes weights and derivative matrix
   and stores the data on a previously allocated jac_quadrature structure.

   \param quad A previously allocated jac_quadrature structure
   \param qtype Quadrature type
   \param a Alpha weight
   \param b Beta weight
   \param ws Workspace 3*quad->Q long used in calculations. If it is null, memory is allocated and at the end released
   \return AN error code or GSL_SUCCESS
   
*/
int 
jac_quadrature_zwd(jac_quadrature *quad, enum jac_quad_type qtype, double a, double b, double *ws)
{
  int allocated = 0;
  if (ws == NULL){
    // Try to allocate the workspace:
    ws = (double *) malloc(3*quad->Q);
    if (ws == NULL)
      {
	GSL_ERROR("Could not allocate workspace memory", GSL_ENOMEM);
      }
    allocated = 1;
  }
  quad->alpha = a;
  quad->beta = b;
  quad->type = qtype;
  int err=0;
  // Calculates the zeros of the quadrature
  switch (quad->type)
    {
    case JAC_GJ:
      err = jac_zeros_gj(quad->x, quad->Q, quad->alpha, quad->beta);
      if (err)
	{
	  GSL_ERROR("Problem calculating the zeros", err);
	}
      err = jac_weights_gj(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the weightd", err);
	}
      err = jac_diffmat_gj(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the differentiation matrix", err);
	}
      break;
    case JAC_GLJ:
      err = jac_zeros_glj(quad->x, quad->Q, quad->alpha, quad->beta);
      if (err)
	{
	  GSL_ERROR("Problem calculating the zeros", err);
	}
      err = jac_weights_glj(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the weightd", err);
	}
      err = jac_diffmat_glj(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the differentiation matrix", err);
	}
      break;
    case JAC_GRJM:
      err = jac_zeros_grjm(quad->x, quad->Q, quad->alpha, quad->beta);
      if (err)
	{
	  GSL_ERROR("Problem calculating the zeros", err);
	}
      err = jac_weights_grjm(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the weightd", err);
	}
      err = jac_diffmat_grjm(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the differentiation matrix", err);
	}
      break;
    case JAC_GRJP:
      err = jac_zeros_grjp(quad->x, quad->Q, quad->alpha, quad->beta);
      if (err)
	{
	  GSL_ERROR("Problem calculating the zeros", err);
	}
      err = jac_weights_grjp(quad->x, quad->w, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the weightd", err);
	}
      err = jac_diffmat_grjp(quad->x, quad->D, quad->Q, quad->alpha, quad->beta, ws);
      if (err)
	{
	  GSL_ERROR("Problem calculating the differentiation matrix", err);
	}
      break;
    default:
      GSL_ERROR("Illegal quadrature type", err);
    }


  if (allocated)
    {
      free(ws);
    }

  return GSL_SUCCESS;
}



/**
   Allocates memory for interpolation matrix for a jac_quadrature structure. 
   It also copies the interpolation points.

   \param quad A previously allocated jac_quadrature structure
   \param npoints Number of interpolation points
   \param xp Interpolation points
   \return Error code or GSL_SUCCESS
 */
int 
jac_interpmat_alloc(jac_quadrature *quad, int npoints, double *xp)
{
  if (npoints < 1)
    {
      GSL_ERROR("The number of interpolating points should be at least 1", GSL_EINVAL);
    }

  // Allocate memory for the interpolation matrix and points
  quad->xp = (double *) malloc( sizeof(double) * (npoints + npoints*quad->Q));
  if (!quad->xp)
    {
      GSL_ERROR("Memory for interpolation matrix could not be allocated", GSL_ENOMEM);
    }

  quad->Imat = quad->xp + npoints;
  int i;
  for (i = 0; i < npoints; ++i)
    quad->xp[i] = xp[i];
  

  
  int err=0;
  
  quad->np = npoints;

  switch (quad->type)
    {
    case JAC_GJ:
      err = jac_interpmat_gj(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case JAC_GLJ:
      err = jac_interpmat_glj(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case JAC_GRJM:
      err = jac_interpmat_grjm(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    case JAC_GRJP:
      err = jac_interpmat_grjp(quad->Imat, quad->xp, npoints, quad->x, quad->Q, quad->alpha, quad->beta);
      break;
    default:
      GSL_ERROR("Illegal quadrature type", err);
    }

  return err;
}


/**
   Frees memory used by jac_ionterpmat_alloc

   \param quad An allocated jac_quadrature structure.
 */
void
jac_interpmat_free(jac_quadrature *quad)
{
  if (quad->xp){
    free(quad->xp);
  }
  quad->np = 0;
  quad->xp = NULL;
  quad->Imat = NULL;
}






/** Aproximates the integral of a function using a a quadrature.

  \f[
   \int_{-1}^1 f(x) dx \approx \sum_{i=0}^{Q-1} w_i f(x_i)
  \f]

  Where \f$x_i\f$ and \f$w_i\f$ are given by the quadrature method emplyed.

  \param quad A structure containing quadrature information
  \param f Value of function at quadrature points
  \return \f$\int_{-1}^{1} f(x) dx \approx \sum_{i=0}^{Q-1} w_i f(x_i)\f$
 */
double 
jac_integrate(jac_quadrature *quad, double *f)
{
    double sum = cblas_ddot(quad->Q, quad->w, 1, f, 1);
    return sum;
}
    

/** Calculates the derivative of a function known at quadrature points given the derivative matrix
    The derivative matrix should have been calculated before using one of the functions jac_**_diffmat.

    The derivative is calculated according to the following equation:
    
    \f[
       \left.\frac{du(x)}{dx}\right|_{x=x_i} = \sum_{j=0}^{Q-1} D_{ij} u(x_j)
    \f]

    In this equation,
    \f[
      D_{ij} = \left.\frac{dh_j(x)}{dx}\right|_{x=x_i}
    \f]

     where \f$h_j(x)\f$ is the Lagrange polynomial through the jth quadrature node
     
 \param quad A structure containing quadrature information
 \param f Value of function at quadrature points
 \param d The Estimated derivative at quadrature points: 
 \return GSL_SUCCESS if everything was ok. Otherwise return an error code
 */
int
jac_differentiate(jac_quadrature *quad, double *f, double *d)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, quad->Q, quad->Q, 1.0, quad->D, quad->Q, f,1, 0.0, d, 1);
    return GSL_SUCCESS;
}


/** Interpolates a function given at quadrature points using the interpolation matrix.
    The interpolation matrix is given by the following equation:
    \f[
      I_{ij} = h_j(x_i)
    \f]

    where \f$x_i\f$ are the points where the function should be interpolated.
    With this matrix, the interpolation consists of the following operation:
    \f[
       f^{interpolated}(x_i) = \sum_{j=0}^{Q-1} f_j h_j(x_i)
    \f]

   \param quad Quadrature information (including interpolation matrix)
   \param f Value of function at quadrature points
   \param fout Interpolated values
   \return GSL_SUCCESS if everything was ok. Otherwise return an error code
   
*/
int
jac_interpolate(jac_quadrature *quad, double *f, double *fout)
{
  if (!quad->np)
    {
      // No interpolation data
      GSL_ERROR("No interpolation info was setup", GSL_EINVAL);
    }
  
  cblas_dgemv(CblasRowMajor, CblasNoTrans, quad->np, quad->Q, 1.0, quad->Imat, quad->Q, f,1, 0.0, fout, 1);
  return GSL_SUCCESS;
}
    
