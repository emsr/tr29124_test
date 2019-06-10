/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <wrap_boost.h>

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 8 + std::cout.precision();

  for (int ih = -500; ih <= +500; ++ih)
    {
      auto h = ih * 0.01;
      std::cout << '\n';
      for (int ia = 0; ia <= +500; ++ia)
	{
	  auto a = ia * 0.01;
	  std::cout << ' ' << std::setw(width) << h
		    << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << __gnu_cxx::owens_t(h, a)
		    << ' ' << std::setw(width) << beast::owens_t(h, a)
		    << ' ' << std::setw(width) << __gnu_cxx::owens_t(h, a) - beast::owens_t(h, a)
		    << '\n';
	}
    }
}
