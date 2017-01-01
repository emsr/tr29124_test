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
#include "gauss_hermite_integrate.h"

#define WKB_TOL (1.0e+04 * TEST_SQRT_TOL0)

template<typename _Tp>
  struct my_f
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp d;
    _Tp e;

    _Tp
    operator()(_Tp x)
    { return (((a * x + b) * x + c ) * x + d ) * x + e; }
  };

template<typename _Tp>
  struct my_fc
  {
    _Tp a;
    _Tp b;
    _Tp c;

    _Tp
    operator()(_Tp x)
    { return a * std::cos(b * x + c); }
  };

int
test_hermite_quad()
{
  my_f<double> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
  my_fc<double> cparams{ 1.0, 1.0, 0.0 };

//  pfoo2 = __gnu_test::gauss_hermite_prob_integrate(params, 3);
//  TEST_SF_VAL(sa, x,  +0.0,   7.01855916896680140676414279747, TEST_TOL0);
//  pfooc = __gnu_test::gauss_hermite_prob_integrate(cparams, 29);
//  TEST_SF_VAL(sa, x,  +0.0,   1.52034690106628080561194014675, TEST_TOL0);

  auto xfoo2 = __gnu_test::gauss_hermite_integrate(params, 3);
  TEST_SF_VAL(sa, x,  +0.0,   2.96886020026673934572443053460, TEST_TOL0);
  xfooc = __gnu_test::gauss_hermite_integrate(cparams, 29);
  TEST_SF_VAL(sa, x,  +0.0,   1.38038844704314297477341524673, TEST_TOL0);

  // printf("status= %d\n", s);
  return s;
}
