/**
 *
 */

#include <lgamma_inv.tcc>

#include <iostream>
#include <cmath>

int
main()
{
  std::cout << std::showpoint;
  for (int i = 16; i <= 100; ++i)
    {
      auto x = 0.1 * i;
      auto y = std::lgamma(x);
      auto xx = lgamma_inv(y);
      std::cout << ' ' << x
		<< ' ' << y
		<< ' ' << xx
		<< ' ' << xx - x
		<< '\n';
    }
}
