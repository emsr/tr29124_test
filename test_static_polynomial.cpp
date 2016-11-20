/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_static_polynomial test_static_polynomial.cpp -lquadmath
./test_static_polynomial
*/

//  Get past a bug....
// $HOME/bin/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__=0 -o test_static_polynomial test_static_polynomial.cpp

#include "static_polynomial.h"
#include <iostream>
#include <complex>
#include <sstream>

int
main()
{
  float aa[5]{1.0f, 2.0f, -1.5f, 0.2f, -0.1f};

  __gnu_cxx::_StaticPolynomial<float, 5> p(aa);
  //__gnu_cxx::_StaticPolynomial<float, 5> a(1.0f, 2.0f, -1.5f, 0.2f, -0.1f);
}
