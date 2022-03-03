/**
 *
 */

#include <cmath>
#include <cfenv>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "sf_elementary.h"

/*
 * Burkhardt driver.
 */
void
driver(std::complex<double> z)
{
  std::complex <double> condn;
  std::complex <double> err;
  std::complex <double> res;
  std::complex <double> res_ult;
  std::complex <double> w;

  std::cout << "\n";
  std::cout << "DRIVER:\n";
  std::cout << "  Demonstrate simple and extended Wright Omega evaluators.\n";

  // Simple evaluator.
  w = emsr::wright_omega(z);

  std::cout << "\n";
  std::cout << "  Calling:\n";
  std::cout << "    w = wrightomega(z);\n";
  std::cout << "  returns:\n";
  std::cout << "    w = omega(" << std::real ( z ) 
       << ", " << std::imag ( z ) 
       << ") =  ( " << std::real ( w )
       << ", " << std::imag ( w ) << ")\n";

  // Extended evaluator.
  w = emsr::detail::wright_omega(z, err, res, condn);

  std::cout << "\n";
  std::cout << "  Calling:\n";
  std::cout << "    wrightomega_ext ( z, w, e, r, condest );\n";
  std::cout << "  returns:\n";
  std::cout << "    w = omega(" << std::real ( z ) 
       << ", " << std::imag ( z ) 
       << ") =  ( " << std::real ( w ) 
       << ", " << std::imag ( w ) << ")\n";
  std::cout << "  e = last update step = ( " << std::real ( err ) 
       << ", " << std::imag ( err ) << ")\n";
  std::cout << "  r = penultimate residual = ( " << std::real ( res ) 
       << ", " << std::imag ( res ) << ")\n";
  std::cout << "  condest = condition number estimate = ( " << std::real ( condn ) 
       << ", " << std::imag ( condn ) << ")\n";

  // Calculate and print ultimate residual.
  res_ult = ( 2.0 * w * w - 8.0 * w - 1.0 ) 
    / std::pow ( 1.0 + w, 6.0 ) * std::pow ( res, 4.0 );
  std::cout << "\n";
  std::cout << "  ultimate residual = ( " << std::real ( res_ult )
       << ", " << std::imag ( res_ult ) << ")\n";

  return;
}

/*
 * Make sure we didn't blow anything.
 */
void
burkhadt_test()
{
  double a;
  double b;
  double pi = M_PI;
  std::complex <double> z;

  //timestamp ( );
  std::cout << "\n";
  std::cout << "TOMS917_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the TOMS917 library.\n";

  a = 0.0;
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0;
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0 + exp ( 1.0 );
  b = 0.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 1.0;
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 1.0;
  b = - pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 2.0 + log ( 2.0 );
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = - 2.0 + log ( 2.0 );
  b = - pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 0.0;
  b = 1.0 + pi / 2.0;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 0.0;
  b = pi;
  z = std::complex<double>(a, b);
  driver ( z );

  a = 1.0;
  b = 1.0;
  z = std::complex<double>(a, b);
  driver ( z );
};

void
test_burkhardt_boundary()
{
  std::complex<double> cond;
  std::complex<double> e;
  double exp_num = 160.0;
  std::string filename = "results.txt";
  std::ofstream fp;
  int i;
  int n = 100;
  double pi = M_PI;
  std::complex<double> r;
  double td;
  std::complex<double> w;
  double x[2];
  double y[2];
  std::complex<double> z;

  std::cout << "\n";
  std::cout << "TEST_BOUNDARY:\n";
  std::cout << "  Test wrightomega_ext() near approximation region boundaries.\n";
  std::cout << "  Store results in a file for comparison with benchmark data.\n";

  fp.open ( filename.c_str ( ) );

  //  We want trailing zeros, to make comparison with benchmark easier.
  fp << std::fixed;

  //  Region 1;
  //  x=(-2.0,1.0] ,y=2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( 2.0 * pi, -1.0 );
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 2;
  //  x=1.0 ,y=(1.0,2*pi)
  x[0] = 1.0;
  y[0] = std::nextafter ( 1.0, 2.0 );
  y[1] = std::nextafter ( 2.0 * pi, 1.0 );
  td = - ( y[1] - y[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 3;
  //  x=(-2.0,1.0] ,y=1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = - ( x[1] - x[0] ) / double(n);
  y[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 4;
  //  x=-2.0 ,y=(1.0,2*pi)
  y[0] = std::nextafter ( 1.0, 2.0 );
  y[1] = std::nextafter ( 2.0 * pi, - 1.0 );
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 5;
  //  x=(-2.0,1.0] ,y=-2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 6;
  //  x=1.0, y=(-2*pi,-1.0)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - 1.0, - 2.0 );
  td = ( y[1] - y[0] ) / double(n);
  x[0] = 1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 7;
  //  x=(-2.0,1.0] ,y=-1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[0] - x[1] ) / double(n);
  y[0] = std::nextafter ( - 1.0, - 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
    }

  //  Region 8;
  //  x=-2.0 ,y=(-2*pi,-1.0)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - 1.0, -2.0 );
  td = ( y[0] - y[1] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 9;
  //  x=-2.0 y=[-1.0,1.0]
  y[0] = - 1.0;
  y[1] = 1.0;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 10;
  //  x=(-2.0,1.0] y=1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[1] - x[0] ) / double(n);
  y[0] = 1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 11
  //  x=1.0 y=[1.0,pi]
  y[0] = 1.0;
  y[1] = pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );
  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 12
  //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //  (on inside)
  td = pi / double(n);
  x[0] = pi / 2.0;

  for ( i = 0; i < n; i++ )
  {
    double a = std::nextafter ( pi, -1.0 ) * cos ( x[0] - td *double(i) ) 
             + std::nextafter ( 1.0, - 1.0 );
    double b = std::nextafter ( pi, -1.0 ) * sin ( x[0] - td *double(i) );
    z = std::complex<double> ( a, b );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 13
  //  x=1.0 y=[-pi,-1.0]
  y[0] = - pi;
  y[1] = - 1.0;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 14
  //  x=(-2.0,1.0] y=-1.0
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  td = ( x[1] - x[0] ) / double(n);
  y[0] = -1.0;

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 15
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( n - 1 - i ) / exp_num ), HUGE_VAL );
    y[0] = std::nextafter ( pi - 0.75 * ( x[0] + 1.0 ), HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 16
  y[0] = 0.75 + pi;
  y[1] = 2.0 * pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 17
  //  x=(-2.0,1.0] ,y=2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = 2.0 * pi;
  td = ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 18
  //  x=1.0 ,y=(pi,2*pi)
  y[0] = std::nextafter ( pi, 6.0 );
  y[1] = std::nextafter ( 2.0 * pi, 1.0 );
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 19
  //  (x-0.1e1)*(x-0.1e1)+y*y=pi*pi)
  //  (on outside)
  td = pi / double( n - 1 );
  y[0] = pi / 2.0;

  for ( i = 0; i < n; i++ )
  {
    y[1] = pi * sin ( y[0] - td * i );
    x[0] = sqrt ( pi * pi - y[1] * y[1] ) + 1.0;
    if ( y[1] < 0 )
    {
      z = std::complex<double> ( std::nextafter ( x[0], HUGE_VAL ), std::nextafter ( y[1], - HUGE_VAL ) );
    }
    else
    {
      z = std::complex<double> ( std::nextafter ( x[0], HUGE_VAL ), std::nextafter ( y[1], HUGE_VAL ) );
    } 
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 20;
  //  x=1.0 ,y=(-2*pi,-pi)
  y[0] = std::nextafter ( - 2.0 * pi, 1.0 );
  y[1] = std::nextafter ( - pi, - 6.0 );
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( 1.0, 2.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 21;
  //  x=(-2.0,1.0] ,y=-2*pi
  x[0] = std::nextafter ( - 2.0, 1.0 );
  x[1] = 1.0;
  y[0] = std::nextafter ( - 2.0 * pi, - 7.0 );
  td = - ( x[1] - x[0] ) / double(n);

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[1] + td *double(i), y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 22
  y[0] = - 0.75 - pi;
  y[1] = - 2.0 * pi;
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 23
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( i ) / exp_num ), HUGE_VAL );
    y[0] = std::nextafter ( - pi + 0.75 * ( x[0] + 1.0 ), - HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 24
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = std::nextafter ( - 1.0 - exp (double( n - 1 - i ) / exp_num ), - HUGE_VAL );
    y[0] = std::nextafter ( - pi + 0.75 * ( x[0] + 1.0 ), HUGE_VAL );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 25
  y[0] = - pi;
  y[1] = - 0.75 - pi;
  td = - ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[1] + td * double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 26
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( - pi, -7.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 27
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( - pi, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( n - 1 - i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 28
  y[0] = std::nextafter ( - pi, 1.0 );
  y[1] = pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, -3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 29
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( pi, 1.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 30
  //  x=(-inf,-2) y=pi^+
  y[0] = std::nextafter ( pi, 7.0 );

  for ( i = 0; i < n; i++ )
  {
    x[0] = - 1.0 - exp (double( n - 1 - i ) / exp_num );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 31
  y[0] = std::nextafter ( pi, 7.0 );
  y[1] = 0.75 + pi;
  td = ( y[1] - y[0] ) / double(n);
  x[0] = std::nextafter ( - 2.0, - 3.0 );

  for ( i = 0; i < n; i++ )
  {
    z = std::complex<double> ( x[0], ( y[0] + td *double(i) ) );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Region 32
  //  x=(-inf,-2) y=pi^+
  for ( i = 0; i < n; i++ )
  {
    x[0] = -1.0 - exp (double( n - 1 - i ) / exp_num );
    y[0] = std::nextafter ( pi - 0.75 * ( x[0] + 1.0 ), 0.1 );
    z = std::complex<double> ( x[0], y[0] );
    w = emsr::detail::wright_omega(z, e, r, cond );
    fp << std::real ( z ) << " " 
       << std::imag ( z ) << " " 
       << std::real ( w ) << " "
       << std::imag ( w ) << "\n"; 
  }

  //  Terminate.
  //
  fp.close ( );

  std::cout << "\n";
  std::cout << "TEST_BOUNDARY:\n";
  std::cout << "  Results saved in file '" << filename << "'\n";

  return;
}

/**
 *
 */
template<typename Tp>
  void
  test_wright_omega(Tp proto = Tp{})
  {
    using Cmplx = std::complex<Tp>;
    const std::string path = "test_wright_omega.txt";
    std::ofstream out(path);

    out.precision(emsr::digits10(proto));
    out << std::showpoint << std::scientific;
    auto w = 8 + out.precision();

    out << '\n';
    for (int i = -100; i <= 100; ++i)
      {
        out << '\n';
	auto x = i * Tp{0.1L};
        for (int j = -100; j <= +100; ++j)
          {
            auto y = j * Tp{0.1L};
	    auto wo = emsr::wright_omega(Cmplx(x, y));
	    out << ' ' << std::setw(w) << x
		<< ' ' << std::setw(w) << y
                << ' ' << std::setw(w) << std::real(wo)
                << ' ' << std::setw(w) << std::imag(wo)
		<< '\n';
          }
      }
  }

/**
 *
 */
template<typename Tp>
  void
  test_lambert_wm(Tp proto = Tp{})
  {
    const std::string path = "test_lambert_wm.txt";
    std::ofstream out(path);

    out.precision(emsr::digits10(proto));
    out << std::showpoint << std::scientific;
    auto w = 8 + out.precision();

    const auto s_1de = emsr::inv_e_v<Tp>;

    out << '\n';
    for (int i = 0; i <= 250; ++i)
      {
	auto x = -s_1de + i * Tp{0.002L};
        if (x >= Tp{0})
          break;
        auto y = emsr::lambert_wm(x);
	out << ' ' << std::setw(w) << x
	    << ' ' << std::setw(w) << y
            << '\n';
      }
  }

/**
 *
 */
template<typename Tp>
  void
  test_lambert_wp(Tp proto = Tp{})
  {
    const std::string path = "test_lambert_wp.txt";
    std::ofstream out(path);

    out.precision(emsr::digits10(proto));
    out << std::showpoint << std::scientific;
    auto w = 8 + out.precision();

    const auto s_1de = emsr::inv_e_v<Tp>;

    out << '\n';
    for (int i = 0; i <= 400; ++i)
      {
	auto x = -s_1de + i * Tp{0.01L};
        auto y = emsr::lambert_wp(x);
	out << ' ' << std::setw(w) << x
	    << ' ' << std::setw(w) << y
            << '\n';
      }
  }

int
main()
{
  burkhadt_test();
  test_burkhardt_boundary();

  test_wright_omega<double>();

  test_lambert_wm<double>();

  test_lambert_wp<double>();
}
