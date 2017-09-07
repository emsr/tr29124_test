/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I.. -o test_gauss_laguerre test_gauss_laguerre.cpp
*/

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gauss_laguerre_integrate.h"

template<typename _Tp>
  struct my_f
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp d;
    _Tp e;

    _Tp
    operator()(_Tp x)
    { return (((a * x + b) * x + c ) * x + d ) * x + e; }
  };

template<typename _Tp>
  struct my_fc
  {
    _Tp a;
    _Tp b;
    _Tp c;

    _Tp
    operator()(_Tp x)
    { return a * std::cos(b * x + c); }
  };

void
test_laguerre_quad()
{
  my_f<double> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
  my_fc<double> cparams{ 1.0, 1.0, 0.0 };

  //pfoo2 = __gnu_ext::gauss_laguerre_prob_integrate<double>(params, 3);
  //pfooc = __gnu_ext::gauss_laguerre_prob_integrate<double>(cparams, 29);

  auto xfoo2 = __gnu_ext::gauss_laguerre_integrate<double>(params, 3, 0.0);
  std::cout << "xfoo2 = " << xfoo2 << '\n';
  std::cout << "delta = " << xfoo2 - 0.0L << '\n';
  auto xfooc = __gnu_ext::gauss_laguerre_integrate<double>(cparams, 29, 0.0);
  std::cout << "xfooc = " << xfooc << '\n';
  std::cout << "delta = " << xfooc - 0.0L << '\n';
}
