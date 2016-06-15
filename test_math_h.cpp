/*
$HOME/bin_tr29124/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_math_h test_math_h.cpp
*/

#include <math.h>

int
main()
{
  int n = 3;
  double x = 4.5;
  double jn = sph_bessel(n, x);
}
