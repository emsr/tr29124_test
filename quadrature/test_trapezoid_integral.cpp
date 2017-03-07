
#include <iostream>
#include <iomanip>
#include <cmath>

#include "trapezoid_integral.h"

int
main()
{
  auto fun = [](long double x) ->long double { return std::sin(x); };
  using fun_t = decltype(fun);
  __gnu_test::trapezoid_integral<fun_t, long double> tq(fun, 0.0l, 3.14159l, 0.0000001l);
  std::cout << tq() << '\n';
}
