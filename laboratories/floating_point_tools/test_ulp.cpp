/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>

template<typename Tp>
  Tp
  ulp(Tp x)
  {
    const int digs = std::numeric_limits<Tp>::digits;
    const int min_exp = std::numeric_limits<Tp>::min_exponent;
    int exp;
    /*Tp frac = */std::frexp(std::abs(x), &exp);
    return ldexp(Tp{1}, std::max(exp, min_exp) - digs + 1);
  }

int
main()
{
  auto a = ulp(-123.457);
  std::cout << "ulp(-123.457) = " << a << '\n';
}
