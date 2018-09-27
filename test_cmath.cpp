/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_cmath test_cmath.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
./test_cmath > test_cmath.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_cmath test_cmath.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
./test_cmath > test_cmath.txt
*/

#include <iostream>
#include <cmath>

int
main()
{
  double x = 0.5;
  double y = -1.0;
  double z = 1.5;
  int i = 0;
  long l = 0;

  std::cout << "  x = " << x  << '\n';
  std::cout << "  y = " << y << '\n';
  std::cout << "  z = " << z << '\n';
  std::cout << "  acos( x ) = " << std::acos( x ) << '\n';
  std::cout << "  acosh( z ) = " << std::acosh( z ) << '\n';
  std::cout << "  asin( x ) = " << std::asin( x ) << '\n';
  std::cout << "  asinh( x ) = " << std::asinh( x ) << '\n';
  std::cout << "  atan( x ) = " << std::atan( x ) << '\n';
  std::cout << "  atan2( x, y ) = " << std::atan2( x, y ) << '\n';
  std::cout << "  atanh( x ) = " << std::atanh( x ) << '\n';
  std::cout << "  cbrt( x ) = " << std::cbrt( x ) << '\n';
  std::cout << "  ceil( x ) = " << std::ceil( x ) << '\n';
  std::cout << "  copysign( x, y ) = " << std::copysign( x, y ) << '\n';
  std::cout << "  cos( x ) = " << std::cos( x ) << '\n';
  std::cout << "  cosh( x ) = " << std::cosh( x ) << '\n';
  std::cout << "  erf( x ) = " << std::erf( x ) << '\n';
  std::cout << "  exp( x ) = " << std::exp( x ) << '\n';
  std::cout << "  exp2( x ) = " << std::exp2( x ) << '\n';
  std::cout << "  expm1( x ) = " << std::expm1( x ) << '\n';
  std::cout << "  fabs( x ) = " << std::fabs( x ) << '\n';
  std::cout << "  fdim( x, y ) = " << std::fdim( x, y ) << '\n';
  std::cout << "  floor( x ) = " << std::floor( x ) << '\n';
  std::cout << "  fma( x, y, z ) = " << std::fma( x, y, z ) << '\n';
  std::cout << "  fmax( x, y ) = " << std::fmax( x, y ) << '\n';
  std::cout << "  fmin( x, y ) = " << std::fmin( x, y ) << '\n';
  std::cout << "  fmod( x, y ) = " << std::fmod( x, y ) << '\n';
  std::cout << "  frexp( x, &i ) = " << std::frexp( x, &i ) << '\n';
  std::cout << "  hypot( x, y ) = " << std::hypot( x, y ) << '\n';
  std::cout << "  ilogb( x ) = " << std::ilogb( x ) << '\n';
  std::cout << "  ldexp( x, i ) = " << std::ldexp( x, i ) << '\n';
  std::cout << "  lgamma( x ) = " << std::lgamma( x ) << '\n';
  std::cout << "  llrint( x ) = " << std::llrint( x ) << '\n';
  std::cout << "  llround( x ) = " << std::llround( x ) << '\n';
  std::cout << "  log( x ) = " << std::log( x ) << '\n';
  std::cout << "  log10( x ) = " << std::log10( x ) << '\n';
  std::cout << "  log2( x ) = " << std::log2( x ) << '\n';
  std::cout << "  logb( x ) = " << std::logb( x ) << '\n';
  std::cout << "  lrint( x ) = " << std::lrint( x ) << '\n';
  std::cout << "  lround( x ) = " << std::lround( x ) << '\n';
  std::cout << "  nearbyint( x ) = " << std::nearbyint( x ) << '\n';
  std::cout << "  nextafter( x, y ) = " << std::nextafter( x, y ) << '\n';
  std::cout << "  nexttoward( x, y ) = " << std::nexttoward( x, y ) << '\n';
  std::cout << "  pow( x, y ) = " << std::pow( x, y ) << '\n';
  std::cout << "  remainder( x, y ) = " << std::remainder( x, y ) << '\n';
  std::cout << "  remquo( x, y, &i ) = " << std::remquo( x, y, &i ) << '\n';
  std::cout << "  rint( x ) = " << std::rint( x ) << '\n';
  std::cout << "  round( x ) = " << std::round( x ) << '\n';
  std::cout << "  scalbln( x, l ) = " << std::scalbln( x, l ) << '\n';
  std::cout << "  scalbn( x, i ) = " << std::scalbn( x, i ) << '\n';
  std::cout << "  sin( x ) = " << std::sin( x ) << '\n';
  std::cout << "  sinh( x ) = " << std::sinh( x ) << '\n';
  std::cout << "  sqrt( x ) = " << std::sqrt( x ) << '\n';
  std::cout << "  tan( x ) = " << std::tan( x ) << '\n';
  std::cout << "  tanh( x ) = " << std::tanh( x ) << '\n';
  std::cout << "  tgamma( x ) = " << std::tgamma( x ) << '\n';
  std::cout << "  trunc( x ) = " << std::trunc( x ) << '\n';

  return 0;
}
