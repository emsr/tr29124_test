
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
//* Hermite polynomials, Hermite functions                               *
//* and their respective arbitrary derivatives                           *
//* Copyright 2011-2013 Konrad Griessinger                               *
//* (konradg(at)gmx.net)                                                 *
//*----------------------------------------------------------------------*

// TODO:
//   - array functions for derivatives of Hermite functions
//   - asymptotic approximation for derivatives of Hermite functions
//   - refine existing asymptotic approximations, especially around x=sqrt(2*n+1) or x=sqrt(2*n+1)*sqrt(2), respectively
//   - asymptotic approximations to prevent overflow also in array functions (?)


#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_math.h>
// #include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

static
double 
gsl_sf_hermite_prob(const int n, const double x)
// Evaluates the probabilists' Hermite polynomial of order n at position x.
// For large n an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used, while for small n the upward recurrence is employed.
{
  // return gsl_sf_hyperg_U(-n/2.,1./2.,x*x/2.)*gsl_sf_pow_int(2,n/2)*(GSL_IS_ODD(n)?M_SQRT2:1);

  if(n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(x == 0.){
    if(GSL_IS_ODD(n)){
      return 0.;
    }
    else{
      double f;
      int j;
      f = (GSL_IS_ODD(n/2)?-1.:1.);
      for(j=1; j < n; j+=2) {
	f*=j;
      }
      return f;
    }
  }
  else if(n == 0) {
    return 1.0;
  }
  else if(n == 1) {
    return x;
  }
/*
  else if(x*x < 4.0*n && n > 100000) {
    // asymptotic formula
    double f = 1.0;
    int j;
    if(GSL_IS_ODD(n)) {
      f=gsl_sf_fact((n-1)/2)*gsl_sf_pow_int(2,n/2)*M_SQRT2/M_SQRTPI;
    }
    else {
      for(j=1; j < n; j+=2) {
	f*=j;
      }
    }
    return f*exp(x*x/4)*cos(x*sqrt(n)-(n%4)*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
    // return f*exp(x*x/4)*cos(x*sqrt(n)-n*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
  }
*/
  else if (n > 10000){
    // Plancherel-Rotach approximation (note: Szego defines the Airy function differently!)
    const double aizero1 = -2.3381074104597670384891972524467; // first zero of the Airy function Ai
    //const double aizero1 = -2.3381074104597670384891972524467354406385401456723878524838544372; // first zero of the Airy function Ai
    double z = fabs(x)*M_SQRT1_2;
    double f = 1.;
    int j;
    for(j=1; j <= n; j++) {
      f*=sqrt(j);
    }
    if (z < sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.)){
      double phi = acos(z/sqrt(2*n+1.));
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(2./n,0.25)/sqrt(M_SQRTPI*sin(phi))*sin(M_PI*0.75+(n/2.+0.25)*(sin(2*phi)-2*phi))*exp(z*z/2.);
    }
    else if (z > sqrt(2*n+1.)-aizero1/pow(8.*n,1/6.)){
      double phi = gsl_acosh(z/sqrt(2*n+1.));
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(0.125/n,0.25)/sqrt(M_SQRTPI*sinh(phi))*exp((n/2.+0.25)*(2*phi-sinh(2*phi)))*exp(z*z/2.);
    }
    else{
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*sqrt(M_SQRTPI)*pow(2.,0.25)*pow(n,-1/12.)*gsl_sf_airy_Ai((z-sqrt(2*n+1.))*pow(8.*n,1/6.),0)*exp(z*z/2.);
    }
  }
  else {
    // upward recurrence: He_{n+1} = x He_n - n He_{n-1} 

    double p_n0 = 1.0;    // He_0(x) 
    double p_n1 = x;    // He_1(x) 
    double p_n = p_n1;
    int j;

    for(j=1; j <= n-1; j++){
      p_n  = x*p_n1-j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
    }

    return p_n;
    }

}

static
double 
gsl_sf_hermite_prob_der(const int m, const int n, const double x)
// Evaluates the m-th derivative of the probabilists' Hermite polynomial of order n at position x.
// The direct formula He^{(m)}_n = n!/(n-m)!*He_{n-m}(x) (where He_j(x) is the j-th probabilists' Hermite polynomial and He^{(m)}_j(x) its m-th derivative) is employed.
{
  if(n < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n < m) {
    return 0.;
  }
  else{
    return gsl_sf_hermite_prob(n-m,x)*gsl_sf_choose(n,m)*gsl_sf_fact(m);
  }
}

static
double
gsl_sf_hermite_phys(const int n, const double x)
// Evaluates the physicists' Hermite polynomial of order n at position x.
// For large n an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used, while for small n the upward recurrence is employed.
{
  // return gsl_sf_hyperg_U(-n/2.,1./2.,x*x)*gsl_sf_pow_int(2,n);

  if(n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(x == 0.){
    if(GSL_IS_ODD(n)){
      return 0.;
    }
    else{
      double f;
      int j;
      f = (GSL_IS_ODD(n/2)?-1.:1.);
      for(j=1; j < n; j+=2) {
	f*=j;
      }
      return gsl_sf_pow_int(2,n/2)*f;
    }
  }
  else if(n == 0) {
    return 1.0;
  }
  else if(n == 1) {
    return 2.0*x;
  }
  /*
  else if(x*x < 2.0*n && n > 100000) {
    // asymptotic formula
    double f = 1.0;
    int j;
    if(GSL_IS_ODD(n)) {
      f=gsl_sf_fact((n-1)/2)*gsl_sf_pow_int(2,n)/M_SQRTPI;
    }
    else {
      for(j=1; j < n; j+=2) {
	f*=j;
      }
      f*=gsl_sf_pow_int(2,n/2);
    }
    return f*exp(x*x/2)*cos(x*sqrt(2.0*n)-(n%4)*M_PI_2)/sqrt(sqrt(1-x*x/2.0/n));
    // return f*exp(x*x/2)*cos(x*sqrt(2.0*n)-n*M_PI_2)/sqrt(sqrt(1-x*x/2.0/n));
  }
  */
  else if (n > 10000){
    // Plancherel-Rotach approximation (note: Szego defines the Airy function differently!)
    const double aizero1 = -2.3381074104597670384891972524467; // first zero of the Airy function Ai
    //const double aizero1 = -2.3381074104597670384891972524467354406385401456723878524838544372; // first zero of the Airy function Ai
    double z = fabs(x);
    double f = 1.;
    int j;
    for(j=1; j <= n; j++) {
      f*=sqrt(j);
    }
    if (z < sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.)){
      double phi = acos(z/sqrt(2*n+1.));
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?M_SQRT2:1.)*gsl_sf_pow_int(2,n/2)*pow(2./n,0.25)/sqrt(M_SQRTPI*sin(phi))*sin(M_PI*0.75+(n/2.+0.25)*(sin(2*phi)-2*phi))*exp(z*z/2.);
    }
    else if (z > sqrt(2*n+1.)-aizero1/pow(8.*n,1/6.)){
      double phi = gsl_acosh(z/sqrt(2*n+1.));
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?M_SQRT2:1.)*gsl_sf_pow_int(2,n/2)*pow(0.125/n,0.25)/sqrt(M_SQRTPI*sinh(phi))*exp((n/2.+0.25)*(2*phi-sinh(2*phi)))*exp(z*z/2.);
    }
    else{
      return f*(GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*(GSL_IS_ODD(n)?M_SQRT2:1.)*sqrt(M_SQRTPI)*gsl_sf_pow_int(2,n/2)*pow(2.,0.25)*pow(n,-1/12.)*gsl_sf_airy_Ai((z-sqrt(2*n+1.))*pow(8.*n,1/6.),0)*exp(z*z/2.);
    }
  }
  else {
    // upward recurrence: H_{n+1} = 2x H_n - 2j H_{n-1} 

    double p_n0 = 1.0;    // H_0(x) 
    double p_n1 = 2.0*x;    // H_1(x) 
    double p_n = p_n1;
    int j;

    for(j=1; j <= n-1; j++){
      p_n  = 2.0*x*p_n1-2.0*j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
    }

    return p_n;
    }

}

static
double 
gsl_sf_hermite_phys_der(const int m, const int n, const double x)
// Evaluates the m-th derivative of the physicists' Hermite polynomial of order n at position x.
// The direct formula H^{(m)}_n = 2**m*n!/(n-m)!*H_{n-m}(x) (where H_j(x) is the j-th physicists' Hermite polynomial and H^{(m)}_j(x) its m-th derivative) is employed.
{
  if(n < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n < m) {
    return 0.;
  }
  return gsl_sf_hermite_phys(n-m,x)*gsl_sf_choose(n,m)*gsl_sf_fact(m)*gsl_sf_pow_int(2,m);
}

static
double
gsl_sf_hermite_func(const int n, const double x)
// Evaluates the Hermite function of order n at position x.
// For large n an approximation depending on the x-range (see Szego, Gabor (1939, 1958, 1967), Orthogonal Polynomials, American Mathematical Society) is used, while for small n the direct formula via the probabilists' Hermite polynomial is applied.
{
  /*
  if (x*x < 2.0*n && n > 100000){
    // asymptotic formula
    double f = 1.0;
    int j;
    // return f*exp(x*x/4)*cos(x*sqrt(n)-n*M_PI_2)/sqrt(sqrt(1-x*x/4.0/n));
    return cos(x*sqrt(2.0*n)-(n%4)*M_PI_2)/sqrt(sqrt(n/M_PI/2.0*(1-x*x/2.0/n)))/M_PI;
  }
  */
  if (x == 0.){
    if (GSL_IS_ODD(n)){
      return 0.;
    }
    else{
      double f;
      int j;
      f = (GSL_IS_ODD(n/2)?-1.:1.);
      for(j=1; j < n; j+=2) {
	f*=sqrt(j/(j+1.));
      }
      return f/sqrt(M_SQRTPI);
    }
  }
  else if (n > 10000){
    // Plancherel-Rotach approximation (note: Szego defines the Airy function differently!)
    const double aizero1 = -2.3381074104597670384891972524467; // first zero of the Airy function Ai
    //const double aizero1 = -2.3381074104597670384891972524467354406385401456723878524838544372; // first zero of the Airy function Ai
    double z = fabs(x);
    if (z < sqrt(2*n+1.)+aizero1/pow(8.*n,1/6.)){
      double phi = acos(z/sqrt(2*n+1.));
      return (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(2./n,0.25)/M_SQRTPI/sqrt(sin(phi))*sin(M_PI*0.75+(n/2.+0.25)*(sin(2*phi)-2*phi));
    }
    else if (z > sqrt(2*n+1.)-aizero1/pow(8.*n,1/6.)){
      double phi = gsl_acosh(z/sqrt(2*n+1.));
      return (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(0.125/n,0.25)/M_SQRTPI/sqrt(sinh(phi))*exp((n/2.+0.25)*(2*phi-sinh(2*phi)));
    }
    else{
      return (GSL_IS_ODD(n)&&(x<0.)?-1.:1.)*pow(2.,0.25)*pow(n,-1/12.)*gsl_sf_airy_Ai((z-sqrt(2*n+1.))*pow(8.*n,1/6.),0);
    }
  }
  else{
    return gsl_sf_hermite_prob(n,M_SQRT2*x)*exp(-x*x/2)/sqrt(M_SQRTPI*gsl_sf_fact(n));
  }
}


int
gsl_sf_hermite_prob_array(const int nmax, const double x, double * result_array)
// Evaluates all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence: He_{n+1} = x He_n - n He_{n-1} 

    double p_n0 = 1.0;    // He_0(x) 
    double p_n1 = x;      // He_1(x) 
    double p_n = p_n1;
    int j;

    result_array[0] = 1.0;
    result_array[1] = x;

    for(j=1; j <= nmax-1; j++){
      p_n  = x*p_n1-j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j+1] = p_n;
    }

    return GSL_SUCCESS;
  }
}



int
gsl_sf_hermite_prob_array_der(const int m, const int nmax, const double x, double * result_array)
// Evaluates the m-th derivative of all probabilists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(nmax < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m == 0) {
    gsl_sf_hermite_prob_array(nmax, x, result_array);
    return GSL_SUCCESS;
  }
  else if(nmax < m) {
    int j;
    for(j=0; j <= nmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(nmax == m) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax] = gsl_sf_fact(m);
    return GSL_SUCCESS;
  }
  else if(nmax == m+1) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax-1] = gsl_sf_fact(m);
    result_array[nmax] = result_array[nmax-1]*(m+1)*x;
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence: He^{(m)}_{n+1} = (n+1)/(n-m+1)*(x He^{(m)}_n - n He^{(m)}_{n-1})

    double p_n0 = gsl_sf_fact(m);    // He^{(m)}_{m}(x) 
    double p_n1 = p_n0*(m+1)*x;      // He^{(m)}_{m+1}(x) 
    double p_n = p_n1;
    int j;

    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }

    result_array[m] = p_n0;
    result_array[m+1] = p_n1;

    for(j=m+1; j <= nmax-1; j++){
      p_n  = (x*p_n1-j*p_n0)*(j+1)/(j-m+1);
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j+1] = p_n;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_hermite_prob_der_array(const int mmax, const int n, const double x, double * result_array)
// Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the probabilists' Hermite polynomial of order n at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(n < 0 || mmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    result_array[0] = 1.0;
    int j;
    for(j=1; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(n == 1 && mmax > 0) {
    result_array[0] = x;
    result_array[1] = 1.0;
    int j;
    for(j=2; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if( mmax == 0) {
    result_array[0] = gsl_sf_hermite_prob(n,x);
    return GSL_SUCCESS;
  }
  else if( mmax == 1) {
    result_array[0] = gsl_sf_hermite_prob(n,x);
    result_array[1] = n*gsl_sf_hermite_prob(n-1,x);
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence

    int k = GSL_MAX_INT(0,n-mmax);
    // Getting a bit lazy here...
    double p_n0 = gsl_sf_hermite_prob(k,x);        // He_k(x) 
    double p_n1 = gsl_sf_hermite_prob(k+1,x);      // He_{k+1}(x) 
    double p_n  = p_n1;
    int j;

    for(j=n+1; j <= mmax; j++){
      result_array[j] = 0.0;
    }

    result_array[GSL_MIN_INT(n,mmax)] = p_n0;
    result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

    for(j=GSL_MIN_INT(mmax,n)-1; j > 0; j--){
      k++;
      p_n  = x*p_n1-k*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j-1] = p_n;
    }

    p_n = 1.0;
    for(j=1; j <= GSL_MIN_INT(n,mmax); j++){
      p_n  = p_n*(n-j+1);
      result_array[j] = p_n*result_array[j];
    }

    return GSL_SUCCESS;
  }
}


static
double
gsl_sf_hermite_prob_series(const int n, const double x, const double * a)
// Evaluates the series sum_{j=0}^n a_j*He_j(x) with He_j being the j-th probabilists' Hermite polynomial.
// For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to probabilists' Hermite polynomials is used.
{
  // CHECK_POINTER(a)

  if(n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    return a[0];
  }
  else if(n == 1) {
    return a[0]+a[1]*x;
  }
  else {
    // downward recurrence: b_n = a_n + x b_{n+1} - (n+1) b_{n+2} 

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;
    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+x*b0-(j+1)*b1;
      b1 = btmp;
    }

    return b0;
  }
}


int
gsl_sf_hermite_phys_array(const int nmax, const double x, double * result_array)
// Evaluates all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = 2.0*x;
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence: H_{n+1} = 2x H_n - 2n H_{n-1} 

    double p_n0 = 1.0;      // H_0(x) 
    double p_n1 = 2.0*x;    // H_1(x) 
    double p_n = p_n1;
    int j;

    result_array[0] = 1.0;
    result_array[1] = 2.0*x;

    for(j=1; j <= nmax-1; j++){
      p_n  = 2.0*x*p_n1-2.0*j*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j+1] = p_n;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_hermite_phys_array_der(const int m, const int nmax, const double x, double * result_array)
// Evaluates the m-th derivative of all physicists' Hermite polynomials up to order nmax at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(nmax < 0 || m < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m == 0) {
    gsl_sf_hermite_phys_array(nmax, x, result_array);
    return GSL_SUCCESS;
  }
  else if(nmax < m) {
    int j;
    for(j=0; j <= nmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(nmax == m) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax] = gsl_sf_pow_int(2,m)*gsl_sf_fact(m);
    return GSL_SUCCESS;
  }
  else if(nmax == m+1) {
    int j;
    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }
    result_array[nmax-1] = gsl_sf_pow_int(2,m)*gsl_sf_fact(m);
    result_array[nmax] = result_array[nmax-1]*2*(m+1)*x;
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence: H^{(m)}_{n+1} = 2(n+1)/(n-m+1)*(x H^{(m)}_n - n H^{(m)}_{n-1})

    double p_n0 = gsl_sf_pow_int(2,m)*gsl_sf_fact(m);    // H^{(m)}_{m}(x) 
    double p_n1 = p_n0*2*(m+1)*x;    // H^{(m)}_{m+1}(x) 
    double p_n = p_n1;
    int j;

    for(j=0; j < m; j++){
      result_array[j] = 0.0;
    }

    result_array[m] = p_n0;
    result_array[m+1] = p_n1;

    for(j=m+1; j <= nmax-1; j++){
      p_n  = (x*p_n1-j*p_n0)*2*(j+1)/(j-m+1);
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j+1] = p_n;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_hermite_phys_der_array(const int mmax, const int n, const double x, double * result_array)
// Evaluates all derivatives (starting from 0) up to the mmax-th derivative of the physicists' Hermite polynomial of order n at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(n < 0 || mmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    result_array[0] = 1.0;
    int j;
    for(j=1; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if(n == 1 && mmax > 0) {
    result_array[0] = 2*x;
    result_array[1] = 1.0;
    int j;
    for(j=2; j <= mmax; j++){
      result_array[j] = 0.0;
    }
    return GSL_SUCCESS;
  }
  else if( mmax == 0) {
    result_array[0] = gsl_sf_hermite_phys(n,x);
    return GSL_SUCCESS;
  }
  else if( mmax == 1) {
    result_array[0] = gsl_sf_hermite_phys(n,x);
    result_array[1] = 2*n*gsl_sf_hermite_phys(n-1,x);
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence

    int k = GSL_MAX_INT(0,n-mmax);
    // Getting a bit lazy here...
    double p_n0 = gsl_sf_hermite_phys(k,x);        // H_k(x) 
    double p_n1 = gsl_sf_hermite_phys(k+1,x);      // H_{k+1}(x) 
    double p_n  = p_n1;
    int j;

    for(j=n+1; j <= mmax; j++){
      result_array[j] = 0.0;
    }

    result_array[GSL_MIN_INT(n,mmax)] = p_n0;
    result_array[GSL_MIN_INT(n,mmax)-1] = p_n1;

    for(j=GSL_MIN_INT(mmax,n)-1; j > 0; j--){
      k++;
      p_n  = 2*x*p_n1-2*k*p_n0;
      p_n0 = p_n1;
      p_n1 = p_n;
      result_array[j-1] = p_n;
    }

    p_n = 1.0;
    for(j=1; j <= GSL_MIN_INT(n,mmax); j++){
      p_n  = p_n*(n-j+1)*2;
      result_array[j] = p_n*result_array[j];
    }

    return GSL_SUCCESS;
  }
}


static
double
gsl_sf_hermite_phys_series(const int n, const double x, const double * a)
// Evaluates the series sum_{j=0}^n a_j*H_j(x) with H_j being the j-th physicists' Hermite polynomial.
// For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to physicists' Hermite polynomials is used.
{
  // CHECK_POINTER(a)

  if(n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    return a[0];
  }
  else if(n == 1) {
    return a[0]+a[1]*2*x;
  }
  else {
    // downward recurrence: b_n = a_n + 2x b_{n+1} - 2(n+1) b_{n+2} 

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;
    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+2*x*b0-2*(j+1)*b1;
      b1 = btmp;
    }

    return b0;
  }
}


int
gsl_sf_hermite_func_array(const int nmax, const double x, double * result_array)
// Evaluates all Hermite functions up to order nmax at position x. The results are stored in result_array.
// Since all polynomial orders are needed, upward recurrence is employed.
{
  // CHECK_POINTER(result_array)

  if(nmax < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(nmax == 0) {
    result_array[0] = exp(-x*x/2.)/sqrt(M_SQRTPI);
    return GSL_SUCCESS;
  }
  else if(nmax == 1) {
    result_array[0] = exp(-x*x/2.)/sqrt(M_SQRTPI);
    result_array[1] = result_array[0]*M_SQRT2*x;
    return GSL_SUCCESS;
  }
  else {
    // upward recurrence: Psi_{n+1} = sqrt(2/(n+1))*x Psi_n - sqrt(n/(n+1)) Psi_{n-1} 

    double p_n0 = exp(-x*x/2.)/sqrt(M_SQRTPI);    // Psi_0(x) 
    double p_n1 = p_n0*M_SQRT2*x;                 // Psi_1(x) 
    double p_n = p_n1;
    int j;

    result_array[0] = p_n0;
    result_array[1] = p_n1;

  for (j=1;j<=nmax-1;j++)
    {
      p_n=(M_SQRT2*x*p_n1-sqrt(j)*p_n0)/sqrt(j+1.);
      p_n0=p_n1;
      p_n1=p_n;
      result_array[j+1] = p_n;
    }

    return GSL_SUCCESS;
  }
}


static
double
gsl_sf_hermite_func_series(const int n, const double x, const double * a)
// Evaluates the series sum_{j=0}^n a_j*Psi_j(x) with Psi_j being the j-th Hermite function.
// For improved numerical stability the Clenshaw algorithm (Clenshaw, C. W. (July 1955). "A note on the summation of Chebyshev series". Mathematical Tables and other Aids to Computation 9 (51): 118–110.) adapted to Hermite functions is used.
{
  // CHECK_POINTER(a)

  if(n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(n == 0) {
    return a[0]*exp(-x*x/2.)/sqrt(M_SQRTPI);
  }
  else if(n == 1) {
    return (a[0]+a[1]*M_SQRT2*x)*exp(-x*x/2.)/sqrt(M_SQRTPI);
  }
  else {
    // downward recurrence: b_n = a_n + sqrt(2/(n+1))*x b_{n+1} - sqrt((n+1)/(n+2)) b_{n+2} 

    double b0 = 0.;
    double b1 = 0.;
    double btmp = 0.;
    int j;

    for(j=n; j >= 0; j--){
      btmp = b0;
      b0  = a[j]+sqrt(2./(j+1))*x*b0-sqrt((j+1.)/(j+2.))*b1;
      b1 = btmp;
    }

    return b0*exp(-x*x/2.)/sqrt(M_SQRTPI);
  }
}


static
double
gsl_sf_hermite_func_der(const int m, const int n, const double x)
// Evaluates the m-th derivative of the Hermite function of order n at position x.
// A summation including upward recurrences is used.
{
  int j;
  double r,b;
  double h0=1.;
  double h1=x;
  double p0=1.;
  double p1=M_SQRT2*x;
  double f=1.;
  // FIXME: asymptotic formula!
  if(m < 0 || n < 0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(m == 0){
    return gsl_sf_hermite_func(n,x);
  }
  else{
  for (j=GSL_MAX_INT(1,n-m+1);j<=n;j++)
    {
      //f*=2.*j;
      f*=sqrt(2.*j);
    }
  //f*=gsl_sf_pow_int(2,GSL_MIN_INT(n,m)/2)*(GSL_IS_ODD(GSL_MIN_INT(n,m))?M_SQRT2:1.);
    //f=sqrt(f);
  if (m>n)
    {
      f=(GSL_IS_ODD(m-n)?-f:f);
      for (j=0;j<GSL_MIN_INT(n,m-n);j++)
	{
	  f*=(m-j)/(j+1.);
	}
    }
  for (j=1;j<=m-n;j++)
    {
      b=x*h1-j*h0;
      h0=h1;
      h1=b;
    }
  b=0.;
  for (j=1;j<=n-m;j++)
    {
      b=(M_SQRT2*x*p1-sqrt(j)*p0)/sqrt(j+1.);
      p0=p1;
      p1=b;
    }
  b=0.;
  r=0.;
  for (j=GSL_MAX_INT(0,m-n);j<=m;j++)
    {
      r+=f*h0*p0;
      b=x*h1-(j+1.)*h0;
      h0=h1;
      h1=b;
      b=(M_SQRT2*x*p1-sqrt(n-m+j+1.)*p0)/sqrt(n-m+j+2.);
      p0=p1;
      p1=b;
      f*=-(m-j)/(j+1.)/sqrt(n-m+j+1.)*M_SQRT1_2;
    }
  return r*exp(-x*x/2.)/sqrt(M_SQRTPI);
}
}
