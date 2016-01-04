
//*----------------------------------------------------------------------*
//* This program is free software: you can redistribute it and/or modify *
//* it under the terms of the GNU General Public License as published by *
//* the Free Software Foundation, either version 3 of the License, or    *
//* (at your option) any later version.                                  *
//*                                                                      *
//* This program is distributed in the hope that it will be useful,      *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
//* GNU General Public License for more details.                         *
//*                                                                      *
//* You should have received a copy of the GNU General Public License    *
//* along with this program. If not, see <http://www.gnu.org/licenses/>. *
//*----------------------------------------------------------------------*

//*----------------------------------------------------------------------*
//* "The purpose of computing is insight, not numbers." - R.W. Hamming   *
//* (Gauss-) Hermite quadrature routines                                 *
//* Copyright 2013-2014 Konrad Griessinger                               *
//* (konradg(at)gmx.net)                                                 *
//*----------------------------------------------------------------------*


#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include "gsl_sf_hermite.h"
// #include <gsl/gsl_sf_hermite.h>

#include "error.h"
#include "eval.h"

static
double 
pow2(int c)
// Small function to calculate integer powers of 2 quickly by bit-shifting when in the standard integer range. Otherwise repeated squaring via gsl_sf_pow_int is used.
{
  if(c<0 && c>-31){
    return 1/((double)(1 << -c));
  }
  else if(c>=0 && c<31){
    return (double)(1 << c);
  }
  else{
  return gsl_sf_pow_int(2,c);
  }
}

double gsl_integration_hermite_prob(const int n, const gsl_function *func)
{
  if(n <= 0) {
    // DOMAIN_ERROR(result);
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 1) {
    // return M_SQRTPI*M_SQRT2*((*func)(0.));
    return M_SQRTPI*M_SQRT2*(GSL_FN_EVAL(func,0.));
  }
  else {
    double s = 0., x = 0.;
    int j = 1;
    // if(GSL_IS_ODD(n)==1) s=((*func)(0.))/gsl_sf_pow_int(gsl_sf_hermite_prob(n-1,0.),2);
    if(GSL_IS_ODD(n)==1) s=(GSL_FN_EVAL(func,0.))/gsl_sf_pow_int(gsl_sf_hermite_prob(n-1,0.),2);
    for (j=1; j<=n/2; j++) {
      x = gsl_sf_hermite_prob_zero(n, j);
      s += ((GSL_FN_EVAL(func,x))+(GSL_FN_EVAL(func,-x)))/gsl_sf_pow_int(gsl_sf_hermite_prob(n-1,x),2);
      // s += (((*func)(x))+((*func)(-x)))/gsl_sf_pow_int(gsl_sf_hermite_prob(n-1,x),2);
    }
    return s*M_SQRTPI*M_SQRT2*(1.-1./n)*gsl_sf_fact(n-2);
  }
}

double gsl_integration_hermite_phys(const int n, const gsl_function *func)
{
  if(n <= 0) {
    // DOMAIN_ERROR(result);
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 1) {
    // return M_SQRTPI*((*func)(0.));
    return M_SQRTPI*(GSL_FN_EVAL(func,0.));
  }
  else {
    double s = 0., x = 0.;
    int j = 1;
    // if(GSL_IS_ODD(n)==1) s=((*func)(0.))/gsl_sf_pow_int(gsl_sf_hermite_prob(n-1,0.),2);
    if(GSL_IS_ODD(n)==1) s=(GSL_FN_EVAL(func,0.))/gsl_sf_pow_int(gsl_sf_hermite_phys(n-1,0.),2);
    for (j=1; j<=n/2; j++) {
      x = gsl_sf_hermite_phys_zero(n, j);
      s += ((GSL_FN_EVAL(func,x))+(GSL_FN_EVAL(func,-x)))/gsl_sf_pow_int(gsl_sf_hermite_phys(n-1,x),2);
      // s += (((*func)(x))+((*func)(-x)))/gsl_sf_pow_int(gsl_sf_hermite_phys(n-1,x),2);
    }
    return s*M_SQRTPI*pow2(n-1)*(1.-1./n)*gsl_sf_fact(n-2);
  }
}
