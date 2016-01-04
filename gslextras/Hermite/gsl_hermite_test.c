
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
//* Simple test function for Hermite polynomials, functions etc.         *
//* Warning: does not output warning or error if deviations occur!       *
//* Copyright 2013 Konrad Griessinger                                    *
//* (konradg(at)gmx.net)                                                 *
//*----------------------------------------------------------------------*


#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_math.h>
// #include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

int testherm(void)
{
  const double aizero1 = -2.3381074104597670384891972524467; // first zero of the Airy function Ai
  double r;
  double x = cos(2*2*M_SQRT2);
  printf("x = %.18e\n", x);
  int j,m,n;

  printf("\n +++ gsl_sf_hermite_prob +++ \n");

  //n = -2;
  //r = gsl_sf_hermite_prob(n,x);
  //printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 1.);

  n = 1;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", x);

  n = 26;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 4.8421e12);

  n = 125;
  r = gsl_sf_hermite_prob(n,0.);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.4e\n", 0.);

  n = 128;
  r = gsl_sf_hermite_prob(n,0.);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.5e\n", 1.64749e107);

  n = 125;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 4.6866e103);

  n = 128;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -1.88503e107);

  n = 10025;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -3.30905e18961);
  printf("reference = -3.30905e18961 (too large)\n");

  n = 10028;
  x = ((sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -7.51548e18967);
  printf("reference = -7.51548e18967 (too large)\n");

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -4.13693e22243);
  printf("reference = -4.13693e22243 (too large)\n");

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2)*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -8.36369e22250);
  printf("reference = -8.36369e22250 (too large)\n");

  n = 10025;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)))*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 7.398864e22273);
  printf("reference = 7.398864e22273 (too large)\n");

  n = 10028;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)))*M_SQRT2;
  r = gsl_sf_hermite_prob(n,x);
  printf("gsl_sf_hermite_prob(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.50713e22281);
  printf("reference = 1.50713e22281 (too large)\n");


  printf("\n +++ gsl_sf_hermite_prob_der +++ \n");
  x = cos(2*2*M_SQRT2);

  //m = -5;
  //n = 128;
  //r =  gsl_sf_hermite_prob_der(m,n,x);
  //printf("gsl_sf_hermite_prob_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  //m = 5;
  //n = -128;
  //r =  gsl_sf_hermite_prob_der(m,n,x);
  //printf("gsl_sf_hermite_prob_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  m = 225;
  n = 128;
  r =  gsl_sf_hermite_prob_der(m,n,x);
  printf("gsl_sf_hermite_prob_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 0.);

  m = 5;
  n = 128;
  r =  gsl_sf_hermite_prob_der(m,n,x);
  printf("gsl_sf_hermite_prob_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -1.43492e112);

  printf("\n +++ gsl_sf_hermite_phys +++ \n");

  //n = -2;
  //r = gsl_sf_hermite_phys(n,x);
  //printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.);

  n = 1;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 2.*x);

  n = 26;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -8.32015e16);

  n = 125;
  r = gsl_sf_hermite_phys(n,0.);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.5e\n", 0.);

  n = 128;
  r = gsl_sf_hermite_phys(n,0.);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.5e\n", 3.03909e126);

  n = 125;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 2.73606e122);

  n = 128;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 3.8616e126);

  n = 10025;
  x = (sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -2.70743e20470);
  printf("reference = -2.70743e20470 (too large)\n");

  n = 10028;
  x = (sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2;
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -1.73922e20477);
  printf("reference = -1.73922e20477 (too large)\n");

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 3.38479e23752);
  printf("reference = 3.38479e23752 (too large)\n");

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.93551e23760);
  printf("reference = 1.93551e23760 (too large)\n");

  n = 10025;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 6.05366e23782);
  printf("reference = 6.05366e23782 (too large)\n");

  n = 10028;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  r = gsl_sf_hermite_phys(n,x);
  printf("gsl_sf_hermite_phys(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 3.48778e23790);
  printf("reference = 3.48778e23790 (too large)\n");


  printf("\n +++ gsl_sf_hermite_phys_der +++ \n");
  x = cos(2*2*M_SQRT2);

  //m = -5;
  //n = 128;
  //r =  gsl_sf_hermite_phys_der(m,n,x);
  //printf("gsl_sf_hermite_phys_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  //m = 5;
  //n = -128;
  //r =  gsl_sf_hermite_phys_der(m,n,x);
  //printf("gsl_sf_hermite_phys_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  m = 225;
  n = 128;
  r =  gsl_sf_hermite_phys_der(m,n,x);
  printf("gsl_sf_hermite_phys_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 0.);

  m = 5;
  n = 128;
  r =  gsl_sf_hermite_phys_der(m,n,x);
  printf("gsl_sf_hermite_phys_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -6.89377e131);

  printf("\n +++ gsl_sf_hermite_func +++ \n");

  //n = -2;
  //r = gsl_sf_hermite_func(n,x);
  //printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 0.54098);

  n = 1;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 0.61984);

  n = 26;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -0.273596);

  n = 125;
  r = gsl_sf_hermite_func(n,0.);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.5e\n", 0.);

  n = 128;
  r = gsl_sf_hermite_func(n,0.);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, 0., r);
  printf("reference = %.5e\n", 0.19928);

  n = 125;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 0.05230);

  n = 128;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 0.18237);

  n = 10025;
  x = (sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -4.20129e-2);

  n = 10028;
  x = (sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.))/2;
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -4.49474e-2);

  n = 10025;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 6.08490e-2);

  n = 10028;
  x = (sqrt(2*n+1.)-(aizero1/pow(8.*n,1/6.))/2);
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 6.08475e-2);

  n = 10025;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.23112e-4);

  n = 10028;
  x = (sqrt(2*n+1.)-2*(aizero1/pow(8.*n,1/6.)));
  r = gsl_sf_hermite_func(n,x);
  printf("gsl_sf_hermite_func(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.23109e-4);


  printf("\n +++ gsl_sf_hermite_func_der +++ \n");
  x = cos(2*2*M_SQRT2);

  //m = -5;
  //n = 128;
  //r =  gsl_sf_hermite_func_der(m,n,x);
  //printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  //m = 5;
  //n = -128;
  //r =  gsl_sf_hermite_func_der(m,n,x);
  //printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  m = 225;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 6.89461e277);

  m = 5;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -8.5143e4);

  x = cos(2*2*M_SQRT2);
  double res[10];

  printf("\n +++ gsl_sf_hermite_prob_array +++ \n");

  //n = -2;
  //gsl_sf_hermite_prob_array(n,x,res);
  //printf("gsl_sf_hermite_prob_array(%d,%g)[%d] = %.18e\n", n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  gsl_sf_hermite_prob_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob(j,x));
  }
  printf("\n");

  n = 1;
  gsl_sf_hermite_prob_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob(j,x));
  }
  printf("\n");

  n = 3;
  gsl_sf_hermite_prob_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob(j,x));
  }
  printf("\n");

  printf("\n +++ gsl_sf_hermite_prob_array_der +++ \n");

  //n = -2;
  //gsl_sf_hermite_prob_array_der(n,x,res);
  //printf("gsl_sf_hermite_prob_array_der(%d,%g)[%d] = %.18e\n", n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 2;
  m = 0;
  gsl_sf_hermite_prob_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(m,j,x));
  }
  printf("\n");

  n = 1;
  m = 2;
  gsl_sf_hermite_prob_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(m,j,x));
  }
  printf("\n");

  n = 2;
  m = 2;
  gsl_sf_hermite_prob_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(m,j,x));
  }
  printf("\n");

  n = 3;
  m = 2;
  gsl_sf_hermite_prob_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(m,j,x));
  }
  printf("\n");

  n = 5;
  m = 2;
  gsl_sf_hermite_prob_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_prob_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(m,j,x));
  }
  printf("\n");


  printf("\n +++ gsl_sf_hermite_prob_der_array +++ \n");

  //n = -2;
  //m = 2;
  //gsl_sf_hermite_prob_der_array(m,n,x,res);
  //printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  m = 2;
  gsl_sf_hermite_prob_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 0;
  gsl_sf_hermite_prob_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 1;
  gsl_sf_hermite_prob_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 3;
  gsl_sf_hermite_prob_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 5;
  gsl_sf_hermite_prob_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_prob_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_prob_der(j,n,x));
  }
  printf("\n");

  printf("\n +++ gsl_sf_hermite_prob_series +++ \n");

  n = 9;
  for(j=0;j<=n;j++){
    res[j] = sqrt(j-aizero1); // recycling the array with arbitrary coefficients
  }

  //n = -2;
  //r = gsl_sf_hermite_prob_series(n,x,res);
  //printf("gsl_sf_hermite_prob_series(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_prob_series(n,x,res);
  printf("gsl_sf_hermite_prob_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.5290871);

  n = 1;
  r = gsl_sf_hermite_prob_series(n,x,res);
  printf("gsl_sf_hermite_prob_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 3.00933);

  n = 7;
  r = gsl_sf_hermite_prob_series(n,x,res);
  printf("gsl_sf_hermite_prob_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", -7.07997e1);

  printf("\n +++ gsl_sf_hermite_phys_array +++ \n");

  //n = -2;
  //gsl_sf_hermite_phys_array(n,x,res);
  //printf("gsl_sf_hermite_phys_array(%d,%g)[%d] = %.18e\n", n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  gsl_sf_hermite_phys_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys(j,x));
  }
  printf("\n");

  n = 1;
  gsl_sf_hermite_phys_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys(j,x));
  }
  printf("\n");

  n = 3;
  gsl_sf_hermite_phys_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys(j,x));
  }
  printf("\n");

  printf("\n +++ gsl_sf_hermite_phys_array_der +++ \n");

  //n = -2;
  //gsl_sf_hermite_phys_array_der(n,x,res);
  //printf("gsl_sf_hermite_phys_array_der(%d,%g)[%d] = %.18e\n", n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 2;
  m = 0;
  gsl_sf_hermite_phys_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(m,j,x));
  }
  printf("\n");

  n = 1;
  m = 2;
  gsl_sf_hermite_phys_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(m,j,x));
  }
  printf("\n");

  n = 2;
  m = 2;
  gsl_sf_hermite_phys_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(m,j,x));
  }
  printf("\n");

  n = 3;
  m = 2;
  gsl_sf_hermite_phys_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(m,j,x));
  }
  printf("\n");

  n = 5;
  m = 2;
  gsl_sf_hermite_phys_array_der(m,n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_phys_array_der(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(m,j,x));
  }
  printf("\n");


  printf("\n +++ gsl_sf_hermite_phys_der_array +++ \n");

  //n = -2;
  //m = 2;
  //gsl_sf_hermite_phys_der_array(m,n,x,res);
  //printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  m = 2;
  gsl_sf_hermite_phys_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 0;
  gsl_sf_hermite_phys_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 1;
  gsl_sf_hermite_phys_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 3;
  gsl_sf_hermite_phys_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(j,n,x));
  }
  printf("\n");

  n = 2;
  m = 5;
  gsl_sf_hermite_phys_der_array(m,n,x,res);
  for(j=0;j<=m;j++){
    printf("gsl_sf_hermite_phys_der_array(%d,%d,%g)[%d] = %.18e\n", m, n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_phys_der(j,n,x));
  }
  printf("\n");

  printf("\n +++ gsl_sf_hermite_phys_series +++ \n");

  n = 9;
  for(j=0;j<=n;j++){
    res[j] = sqrt(j-aizero1); // recycling the array with arbitrary coefficients
  }

  //n = -2;
  //r = gsl_sf_hermite_phys_series(n,x,res);
  //printf("gsl_sf_hermite_phys_series(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_phys_series(n,x,res);
  printf("gsl_sf_hermite_phys_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.5290871);

  n = 1;
  r = gsl_sf_hermite_phys_series(n,x,res);
  printf("gsl_sf_hermite_phys_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 4.48958);

  n = 7;
  r = gsl_sf_hermite_phys_series(n,x,res);
  printf("gsl_sf_hermite_phys_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 4.5477e2);

  printf("\n +++ gsl_sf_hermite_func_array +++ \n");

  //n = -2;
  //gsl_sf_hermite_func_array(n,x,res);
  //printf("gsl_sf_hermite_func_array(%d,%g)[%d] = %.18e\n", n, x, 0, res[0]);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  gsl_sf_hermite_func_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_func_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_func(j,x));
  }
  printf("\n");

  n = 1;
  gsl_sf_hermite_func_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_func_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_func(j,x));
  }
  printf("\n");

  n = 3;
  gsl_sf_hermite_func_array(n,x,res);
  for(j=0;j<=n;j++){
    printf("gsl_sf_hermite_func_array(%d,%g)[%d] = %.18e\n", n, x, j, res[j]);
    printf("reference = %.4e\n", gsl_sf_hermite_func(j,x));
  }
  printf("\n");

  printf("\n +++ gsl_sf_hermite_func_series +++ \n");

  n = 9;
  for(j=0;j<=n;j++){
    res[j] = sqrt(j-aizero1); // recycling the array with arbitrary coefficients
  }

  //n = -2;
  //r = gsl_sf_hermite_func_series(n,x,res);
  //printf("gsl_sf_hermite_func_series(%d,%g) = %.18e\n", n, x, r);
  //printf("reference = ERROR / undefined \n");

  n = 0;
  r = gsl_sf_hermite_func_series(n,x,res);
  printf("gsl_sf_hermite_func_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 0.827199);

  n = 1;
  r = gsl_sf_hermite_func_series(n,x,res);
  printf("gsl_sf_hermite_func_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.5e\n", 1.95967);

  n = 7;
  r = gsl_sf_hermite_func_series(n,x,res);
  printf("gsl_sf_hermite_func_series(%d,%g) = %.18e\n", n, x, r);
  printf("reference = %.4e\n", 2.06062);

  printf("\n +++ gsl_sf_hermite_func_der +++ \n");
  x = cos(2*2*M_SQRT2);

  //m = -5;
  //n = 128;
  //r =  gsl_sf_hermite_func_der(m,n,x);
  //printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  //m = 5;
  //n = -128;
  //r =  gsl_sf_hermite_func_der(m,n,x);
  //printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  //printf("reference = ERROR / undefined \n");

  m = 0;
  n = 127;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -0.0714512);

  m = 1;
  n = 127;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -2.97576);

  m = 2;
  n = 127;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 18.1732);

  m = 225;
  n = 127;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 1.50931e278);

  m = 0;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 0.182367);

  m = 1;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -1.29097);

  m = 2;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", -46.7485);

  m = 225;
  n = 128;
  r =  gsl_sf_hermite_func_der(m,n,x);
  printf("gsl_sf_hermite_func_der(%d,%d,%g) = %.18e\n", m, n, x, r);
  printf("reference = %.5e\n", 6.89461e277);

  return 0;
}
