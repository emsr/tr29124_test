/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_cmath test_cmath.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
./test_cmath > test_cmath.txt
*/

#include <iostream>
#include <tr1/cmath>

#if ! _GLIBCXX_USE_C99_MATH_TR1
#error "Unable to compile because all these functions will not be in std::tr1::"
#endif

int
main()
{
  double x = 0.5;
  double y = -1.0;
  double z = 1.5;
  int i = 0;
  long l = 0;

  std::cout << "  x = " << x  << std::endl;
  std::cout << "  y = " << y << std::endl;
  std::cout << "  z = " << z << std::endl;
  std::cout << "  acos( x ) = " << std::tr1::acos( x ) << std::endl;
  std::cout << "  acosh( z ) = " << std::tr1::acosh( z ) << std::endl;
  std::cout << "  asin( x ) = " << std::tr1::asin( x ) << std::endl;
  std::cout << "  asinh( x ) = " << std::tr1::asinh( x ) << std::endl;
  std::cout << "  atan( x ) = " << std::tr1::atan( x ) << std::endl;
  std::cout << "  atan2( x, y ) = " << std::tr1::atan2( x, y ) << std::endl;
  std::cout << "  atanh( x ) = " << std::tr1::atanh( x ) << std::endl;
  std::cout << "  cbrt( x ) = " << std::tr1::cbrt( x ) << std::endl;
  std::cout << "  ceil( x ) = " << std::tr1::ceil( x ) << std::endl;
  std::cout << "  copysign( x, y ) = " << std::tr1::copysign( x, y ) << std::endl;
  std::cout << "  cos( x ) = " << std::tr1::cos( x ) << std::endl;
  std::cout << "  cosh( x ) = " << std::tr1::cosh( x ) << std::endl;
  std::cout << "  erf( x ) = " << std::tr1::erf( x ) << std::endl;
  std::cout << "  exp( x ) = " << std::tr1::exp( x ) << std::endl;
  std::cout << "  exp2( x ) = " << std::tr1::exp2( x ) << std::endl;
  std::cout << "  expm1( x ) = " << std::tr1::expm1( x ) << std::endl;
  std::cout << "  fabs( x ) = " << std::tr1::fabs( x ) << std::endl;
  std::cout << "  fdim( x, y ) = " << std::tr1::fdim( x, y ) << std::endl;
  std::cout << "  floor( x ) = " << std::tr1::floor( x ) << std::endl;
  std::cout << "  fma( x, y, z ) = " << std::tr1::fma( x, y, z ) << std::endl;
  std::cout << "  fmax( x, y ) = " << std::tr1::fmax( x, y ) << std::endl;
  std::cout << "  fmin( x, y ) = " << std::tr1::fmin( x, y ) << std::endl;
  std::cout << "  fmod( x, y ) = " << std::tr1::fmod( x, y ) << std::endl;
  std::cout << "  frexp( x, &i ) = " << std::tr1::frexp( x, &i ) << std::endl;
  std::cout << "  hypot( x, y ) = " << std::tr1::hypot( x, y ) << std::endl;
  std::cout << "  ilogb( x ) = " << std::tr1::ilogb( x ) << std::endl;
  std::cout << "  ldexp( x, i ) = " << std::tr1::ldexp( x, i ) << std::endl;
  std::cout << "  lgamma( x ) = " << std::tr1::lgamma( x ) << std::endl;
  std::cout << "  llrint( x ) = " << std::tr1::llrint( x ) << std::endl;
  std::cout << "  llround( x ) = " << std::tr1::llround( x ) << std::endl;
  std::cout << "  log( x ) = " << std::tr1::log( x ) << std::endl;
  std::cout << "  log10( x ) = " << std::tr1::log10( x ) << std::endl;
  std::cout << "  log2( x ) = " << std::tr1::log2( x ) << std::endl;
  std::cout << "  logb( x ) = " << std::tr1::logb( x ) << std::endl;
  std::cout << "  lrint( x ) = " << std::tr1::lrint( x ) << std::endl;
  std::cout << "  lround( x ) = " << std::tr1::lround( x ) << std::endl;
  std::cout << "  nearbyint( x ) = " << std::tr1::nearbyint( x ) << std::endl;
  std::cout << "  nextafter( x, y ) = " << std::tr1::nextafter( x, y ) << std::endl;
  std::cout << "  nexttoward( x, y ) = " << std::tr1::nexttoward( x, y ) << std::endl;
  std::cout << "  pow( x, y ) = " << std::tr1::pow( x, y ) << std::endl;
  std::cout << "  remainder( x, y ) = " << std::tr1::remainder( x, y ) << std::endl;
  std::cout << "  remquo( x, y, &i ) = " << std::tr1::remquo( x, y, &i ) << std::endl;
  std::cout << "  rint( x ) = " << std::tr1::rint( x ) << std::endl;
  std::cout << "  round( x ) = " << std::tr1::round( x ) << std::endl;
  std::cout << "  scalbln( x, l ) = " << std::tr1::scalbln( x, l ) << std::endl;
  std::cout << "  scalbn( x, i ) = " << std::tr1::scalbn( x, i ) << std::endl;
  std::cout << "  sin( x ) = " << std::tr1::sin( x ) << std::endl;
  std::cout << "  sinh( x ) = " << std::tr1::sinh( x ) << std::endl;
  std::cout << "  sqrt( x ) = " << std::tr1::sqrt( x ) << std::endl;
  std::cout << "  tan( x ) = " << std::tr1::tan( x ) << std::endl;
  std::cout << "  tanh( x ) = " << std::tr1::tanh( x ) << std::endl;
  std::cout << "  tgamma( x ) = " << std::tr1::tgamma( x ) << std::endl;
  std::cout << "  trunc( x ) = " << std::tr1::trunc( x ) << std::endl;

  return 0;
}
