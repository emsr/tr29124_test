/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_inv_lgamma test_inv_lgamma.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_inv_lgamma > test_inv_lgamma.txt
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
