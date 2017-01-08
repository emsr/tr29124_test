/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_float128 test_float128.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_float128 > test_float128.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_float128 test_float128.cpp -lquadmath
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm> // For clamp

#include <bits/float128_io.h>

int
main()
{
  std::cout.precision(std::numeric_limits<__float128>::max_digits10);
  //auto width = 8 + std::cout.precision();

  auto x = 0.00004472229441850588228136889483397204368247Q;
  auto y = -1.00200035757357000394547575997277994404042345Q;
  auto z = 2.12005432105145066363366054636636034336318985Q;
  int n = 5;
  int exp = 3;
  __float128 m = 4.0Q;
  __float128 b = 2.5Q;
  __float128 xsin, xcos, i;

  //std::cout << std::hexfloat << std::uppercase;

  std::cout << "x               = " << x << '\n';
  std::cout << "y               = " << y << '\n';
  std::cout << "n               = " << n << '\n';
  std::cout << "m               = " << m << '\n';
  std::cout << "b               = " << b << '\n';
  std::cout << "exp             = " << exp << '\n';
  std::cout << "abs(x)          = " << std::abs(x) << '\n';
  std::cout << "acos(x)         = " << std::acos(x) << '\n';
  std::cout << "acosh(x)        = " << std::acosh(x) << '\n';
  std::cout << "asin(x)         = " << std::asin(x) << '\n';
  std::cout << "asinh(x)        = " << std::asinh(x) << '\n';
  std::cout << "atan(x)         = " << std::atan(x) << '\n';
  std::cout << "atanh(x)        = " << std::atanh(x) << '\n';
  std::cout << "atan2(x)        = " << std::atan2(y, x) << '\n';
  std::cout << "cbrt(x)         = " << std::cbrt(x) << '\n';
  std::cout << "ceil(x)         = " << std::ceil(x) << '\n';
  std::cout << "clamp(z, y, x)  = " << std::clamp(z, y, x) << '\n';
  std::cout << "copysign(x, y)  = " << std::copysign(x, y) << '\n';
  std::cout << "cos(x)          = " << std::cos(x) << '\n';
  std::cout << "cosh(x)         = " << std::cosh(x) << '\n';
  std::cout << "exp(x)          = " << std::exp(x) << '\n';
  std::cout << "erf(x)          = " << std::erf(x) << '\n';
  std::cout << "erfc(x)         = " << std::erfc(x) << '\n';
  std::cout << "expm1(x)        = " << std::expm1(x) << '\n';
  std::cout << "fabs(x)         = " << std::fabs(x) << '\n';
  std::cout << "fdim(x)         = " << std::fdim(x, y) << '\n';
  std::cout << "floor(x)        = " << std::floor(x) << '\n';
  std::cout << "fma(m, x, b)    = " << std::fma(m, x, b) << '\n';
  std::cout << "fmax(x, y)      = " << std::fmax(x, y) << '\n';
  std::cout << "fmin(x, y)      = " << std::fmin(x, y) << '\n';
  std::cout << "fmod(x, y)      = " << std::fmod(x, y) << '\n';
  std::cout << "frexp(x, &exp)  = " << std::frexp(x, &exp) << ", exp = " << exp << '\n';
  std::cout << "hypot(x, y)     = " << std::hypot(x, y) << '\n';
  std::cout << "hypot(x, y, z)  = " << std::hypot(x, y, z) << '\n';
  std::cout << "isinf(x)        = " << std::isinf(x) << '\n';
  std::cout << "ilogb(x)        = " << std::ilogb(x) << '\n';
  std::cout << "isnan(x)        = " << std::isnan(x) << '\n';
  std::cout << "j0(x)           = " << std::j0(x) << '\n';
  std::cout << "j1(x)           = " << std::j1(x) << '\n';
  std::cout << "jn(x)           = " << std::jn(n, x) << '\n';
  std::cout << "ldexp(x)        = " << std::ldexp(x, exp) << '\n';
  std::cout << "lgamma(x)       = " << std::lgamma(x) << '\n';
  std::cout << "llrint(x)       = " << std::llrint(x) << '\n';
  std::cout << "llround(x)      = " << std::llround(x) << '\n';
  std::cout << "logb(x)         = " << std::logb(x) << '\n';
  std::cout << "log(x)          = " << std::log(x) << '\n';
  std::cout << "log10(x)        = " << std::log10(x) << '\n';
  std::cout << "log2(x)         = " << std::log2(x) << '\n';
  std::cout << "log1p(x)        = " << std::log1p(x) << '\n';
  std::cout << "lrint(x)        = " << std::lrint(x) << '\n';
  std::cout << "lround(x)       = " << std::lround(x) << '\n';
  std::cout << "modf(x, &i)     = " << std::modf(x, &i) << ", i = " << i << '\n';
  std::cout << "nanq(\"Hi!!!\")   = " << std::nanq("Hi!!!") << '\n';
  std::cout << "nearbyint(x)    = " << std::nearbyint(x) << '\n';
  std::cout << "nextafter(x, y) = " << std::nextafter(x, y) << '\n';
  //std::cout << "powi(x, n)      = " << std::powi(x, n) << '\n';
  std::cout << "remainder(x, y) = " << std::remainder(x, y) << '\n';
  std::cout << "remquo(x, y, n) = " << std::remquo(x, y, &n) << '\n';
  std::cout << "rint(x) 	= " << std::rint(x) << '\n';
  std::cout << "round(x)	= " << std::round(x) << '\n';
  std::cout << "scalbln(x, n)	= " << std::scalbln(x, n) << '\n';
  std::cout << "scalbn(x, n)	= " << std::scalbn(x, n) << '\n';
  std::cout << "signbit(x)	= " << std::signbit(x) << '\n';
  std::sincos(x, &xsin, &xcos);
  std::cout << "sincos(x, &sin, &cos): " << xsin << ", " << xcos << '\n';
  std::cout << "sin(x)		= " << std::sin(x) << '\n';
  std::cout << "sinh(x) 	= " << std::sinh(x) << '\n';
  std::cout << "sqrt(x) 	= " << std::sqrt(x) << '\n';
  std::cout << "tan(x)		= " << std::tan(x) << '\n';
  std::cout << "tanh(x) 	= " << std::tanh(x) << '\n';
  std::cout << "tgamma(x)	= " << std::tgamma(x) << '\n';
  std::cout << "trunc(x)	= " << std::trunc(x) << '\n';
  std::cout << "y0(x)		= " << std::y0(x) << '\n';
  std::cout << "y1(x)		= " << std::y1(x) << '\n';
  std::cout << "yn(n, x)	= " << std::yn(n, x) << '\n';
}
