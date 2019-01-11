/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_tr1_cmath test_tr1_cmath.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
./test_tr1_cmath > test_tr1_cmath.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_tr1_cmath test_tr1_cmath.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
./test_tr1_cmath > test_tr1_cmath.txt
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

  std::cout << "  x = " << x  << '\n';
  std::cout << "  y = " << y << '\n';
  std::cout << "  z = " << z << '\n';
  std::cout << "  acos( x ) = " << std::tr1::acos( x ) << '\n';
  std::cout << "  acosh( z ) = " << std::tr1::acosh( z ) << '\n';
  std::cout << "  asin( x ) = " << std::tr1::asin( x ) << '\n';
  std::cout << "  asinh( x ) = " << std::tr1::asinh( x ) << '\n';
  std::cout << "  atan( x ) = " << std::tr1::atan( x ) << '\n';
  std::cout << "  atan2( x, y ) = " << std::tr1::atan2( x, y ) << '\n';
  std::cout << "  atanh( x ) = " << std::tr1::atanh( x ) << '\n';
  std::cout << "  cbrt( x ) = " << std::tr1::cbrt( x ) << '\n';
  std::cout << "  ceil( x ) = " << std::tr1::ceil( x ) << '\n';
  std::cout << "  copysign( x, y ) = " << std::tr1::copysign( x, y ) << '\n';
  std::cout << "  cos( x ) = " << std::tr1::cos( x ) << '\n';
  std::cout << "  cosh( x ) = " << std::tr1::cosh( x ) << '\n';
  std::cout << "  erf( x ) = " << std::tr1::erf( x ) << '\n';
  std::cout << "  exp( x ) = " << std::tr1::exp( x ) << '\n';
  std::cout << "  exp2( x ) = " << std::tr1::exp2( x ) << '\n';
  std::cout << "  expm1( x ) = " << std::tr1::expm1( x ) << '\n';
  std::cout << "  fabs( x ) = " << std::tr1::fabs( x ) << '\n';
  std::cout << "  fdim( x, y ) = " << std::tr1::fdim( x, y ) << '\n';
  std::cout << "  floor( x ) = " << std::tr1::floor( x ) << '\n';
  std::cout << "  fma( x, y, z ) = " << std::tr1::fma( x, y, z ) << '\n';
  std::cout << "  fmax( x, y ) = " << std::tr1::fmax( x, y ) << '\n';
  std::cout << "  fmin( x, y ) = " << std::tr1::fmin( x, y ) << '\n';
  std::cout << "  fmod( x, y ) = " << std::tr1::fmod( x, y ) << '\n';
  std::cout << "  frexp( x, &i ) = " << std::tr1::frexp( x, &i ) << '\n';
  std::cout << "  hypot( x, y ) = " << std::tr1::hypot( x, y ) << '\n';
  std::cout << "  ilogb( x ) = " << std::tr1::ilogb( x ) << '\n';
  std::cout << "  ldexp( x, i ) = " << std::tr1::ldexp( x, i ) << '\n';
  std::cout << "  lgamma( x ) = " << std::tr1::lgamma( x ) << '\n';
  std::cout << "  llrint( x ) = " << std::tr1::llrint( x ) << '\n';
  std::cout << "  llround( x ) = " << std::tr1::llround( x ) << '\n';
  std::cout << "  log( x ) = " << std::tr1::log( x ) << '\n';
  std::cout << "  log10( x ) = " << std::tr1::log10( x ) << '\n';
  std::cout << "  log2( x ) = " << std::tr1::log2( x ) << '\n';
  std::cout << "  logb( x ) = " << std::tr1::logb( x ) << '\n';
  std::cout << "  lrint( x ) = " << std::tr1::lrint( x ) << '\n';
  std::cout << "  lround( x ) = " << std::tr1::lround( x ) << '\n';
  std::cout << "  nearbyint( x ) = " << std::tr1::nearbyint( x ) << '\n';
  std::cout << "  nextafter( x, y ) = " << std::tr1::nextafter( x, y ) << '\n';
  std::cout << "  nexttoward( x, y ) = " << std::tr1::nexttoward( x, y ) << '\n';
  std::cout << "  pow( x, y ) = " << std::tr1::pow( x, y ) << '\n';
  std::cout << "  remainder( x, y ) = " << std::tr1::remainder( x, y ) << '\n';
  std::cout << "  remquo( x, y, &i ) = " << std::tr1::remquo( x, y, &i ) << '\n';
  std::cout << "  rint( x ) = " << std::tr1::rint( x ) << '\n';
  std::cout << "  round( x ) = " << std::tr1::round( x ) << '\n';
  std::cout << "  scalbln( x, l ) = " << std::tr1::scalbln( x, l ) << '\n';
  std::cout << "  scalbn( x, i ) = " << std::tr1::scalbn( x, i ) << '\n';
  std::cout << "  sin( x ) = " << std::tr1::sin( x ) << '\n';
  std::cout << "  sinh( x ) = " << std::tr1::sinh( x ) << '\n';
  std::cout << "  sqrt( x ) = " << std::tr1::sqrt( x ) << '\n';
  std::cout << "  tan( x ) = " << std::tr1::tan( x ) << '\n';
  std::cout << "  tanh( x ) = " << std::tr1::tanh( x ) << '\n';
  std::cout << "  tgamma( x ) = " << std::tr1::tgamma( x ) << '\n';
  std::cout << "  trunc( x ) = " << std::tr1::trunc( x ) << '\n';

  return 0;
}
