// test_hermite_quad.c

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
//* Test function for Hermite quadrature                                 *
//* Copyright 2013-2014 Konrad Griessinger                               *
//* (konradg(at)gmx.net)                                                 *
//*----------------------------------------------------------------------*


// #include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"
// #include "gsl_sf_hermite.h"
#include "gsl_integration_hermite.h"

// #define PRINT(n) printf("%22.18g  %22.18g  %22.18g  %22.18g\n", F[n], Fp[n], G[n], Gp[n])

#define WKB_TOL (1.0e+04 * TEST_SQRT_TOL0)

struct my_f_params { double a; double b; double c; double d; double e; };

double
my_f (double x, void * p) {
   struct my_f_params * params 
     = (struct my_f_params *)p;
   double a = (params->a);
   double b = (params->b);
   double c = (params->c);
   double d = (params->d);
   double e = (params->e);

   return  ( ( (a * x + b) * x + c ) * x + d ) * x + e;
}

struct my_f_params params = { 0.5, 0., 0., -1., 1.3 };

struct my_fc_params { double a; double b; double c; };

double
my_fc (double x, void * p) {
   struct my_fc_params * cparams 
     = (struct my_fc_params *)p;
   double a = (cparams->a);
   double b = (cparams->b);
   double c = (cparams->c);

   return a*cos(b*x+c);
}

struct my_fc_params cparams = { 1., 1., 0. };

int test_hermite_quad(void)
{
  int s = 0, sa = 0;
  double x;

  gsl_function foo2;
  foo2.function = &my_f;
  foo2.params = &params;

  gsl_function fooc;
  fooc.function = &my_fc;
  fooc.params = &cparams;

  sa = 0;
  x = gsl_integration_hermite_prob(3, &foo2);
  TEST_SF_VAL(sa, x,  +0.0,   7.01855916896680140676414279747, TEST_TOL0);
  x = gsl_integration_hermite_prob(29, &fooc);
  TEST_SF_VAL(sa, x,  +0.0,   1.52034690106628080561194014675, TEST_TOL0);
  gsl_test(sa, "gsl_integration_hermite_prob(n, &foo)");
  s += sa;

  sa = 0;
  x = gsl_integration_hermite_phys(3, &foo2);
  TEST_SF_VAL(sa, x,  +0.0,   2.96886020026673934572443053460, TEST_TOL0);
  x = gsl_integration_hermite_phys(29, &fooc);
  TEST_SF_VAL(sa, x,  +0.0,   1.38038844704314297477341524673, TEST_TOL0);
  gsl_test(sa, "gsl_integration_hermite_phys(n, &foo)");
  s += sa;

  // printf("status= %d\n", s);
  return s;
}
