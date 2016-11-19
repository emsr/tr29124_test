/*
$HOME/bin_tr29124/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -Wno-psabi -I. -o test_math_h test_math_h.cpp -lquadmath
./test_math_h
*/

#include <math.h>

int
main()
{
  unsigned int n = 3;
  double x = 4.5;
  double jn = sph_bessel(n, x);
}
