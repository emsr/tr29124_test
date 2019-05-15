/*
g++ -std=c++14 -Iinclude -o test_ulp test_ulp.cpp
./test_ulp 
*/

#include <ext/ulp.h>
#include <iostream>

int
main()
{
  std::cout << '\n';
  double x = 0.00001;
  for (int i = -5; i <= 5; ++i)
    {
      std::cout << ' ' << i << ' ' << x << ' ' << std::ilogb(x) << '\n';
      x *= 10.0;
    }

  std::cout << '\n';
  x = 0.00001;
  for (int i = -5; i <= 5; ++i)
    {
      std::cout << ' ' << i << ' ' << x << ' ' << __gnu_cxx::ulp(x) << '\n';
      x *= 10.0;
    }
}
