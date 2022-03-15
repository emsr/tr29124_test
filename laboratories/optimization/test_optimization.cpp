
/*
$HOME/bin/bin/g++ -std=c++20 -g -I. -o test_optimization test_optimization.cpp
*/

#include <cmath>
#include <iostream>

#include <optimization.h>

int
main()
{
  auto func1 = [](double x)->double{ return std::cos(x); };
  auto deriv1 = [](double x)->double{ return -std::sin(x); };

  emsr::MinBracket<double> mb(func1, 1.0, 1.5);
  emsr::mean_bracket(func1, mb);
  std::cout << "a.arg = " << mb.a.arg << "; a.val = " << mb.a.val << '\n';
  std::cout << "b.arg = " << mb.b.arg << "; b.val = " << mb.b.val << '\n';
  std::cout << "c.arg = " << mb.c.arg << "; c.val = " << mb.c.val << '\n';

  emsr::MinBracket<double> mbg(func1, 2.0, 3.0, 4.0);
  auto gmin = emsr::golden(func1, mbg, 1.0e-8);
  std::cout << "gmin.arg = " << gmin.arg
            << "; gmin.val = " << gmin.val << '\n';

  auto bmin = emsr::brent(func1, 2.0, 3.0, 4.0, 1.0e-8);
  std::cout << "bmin.arg = " << bmin.arg
            << "; bmin.val = " << bmin.val << '\n';

  auto bdmin = emsr::brent(func1, deriv1, 2.0, 3.0, 4.0, 1.0e-8);
  std::cout << "bdmin.arg = " << bdmin.arg
            << "; bdmin.val = " << bdmin.val
            << "; bdmin.der = " << bdmin.der << '\n';
}
