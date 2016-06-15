// $HOME/bin_tr29124/bin/g++ -std=gnu++1z -o test_pochhammer_l test_pochhammer_l.cpp boost_wrap.cpp
/*
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_pochhammer_l > test_pochhammer_l.txt
*/

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include "boost_wrap.h"

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 8 + std::cout.precision();

  std::vector<double> xv{0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};

  for (int ia = 0; ia <= +500; ++ia)
    {
      auto a = ia * 0.01;
      std::cout << '\n';
      for (int ix = 0; ix < xv.size(); ++ix)
	{
	  auto x = xv[ix];
	  std::cout << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << x
		    << ' ' << std::setw(width) << __gnu_cxx::pochhammer_l(a, x)
		    << ' ' << std::setw(width) << beast::pochhammer_l(a, x)
		    << ' ' << std::setw(width) << __gnu_cxx::pochhammer_l(a, x) - beast::pochhammer_l(a, x)
		    << '\n';
	}
    }
}
