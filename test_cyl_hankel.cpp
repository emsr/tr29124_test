// $HOME/bin_specfun/bin/g++ -g -o test_cyl_hankel test_cyl_hankel.cpp -lquadmath boost_wrap.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_cyl_hankel > test_cyl_hankel.txt

// g++ -std=c++14 -g -o test_cyl_hankel test_cyl_hankel.cpp -lquadmath

// ./test_cyl_hankel > test_cyl_hankel.txt

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include "boost_wrap.h"

template<typename _Tp>
  void
  RunCylHankel1()
  {
    std::cout.precision(std::numeric_limits<double>::digits10);
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    for (int inu = 0; inu <= +500; ++inu)
      {
	auto nu = inu * 0.1;
	std::cout << '\n';
	for (int iz = 1; iz <= +1000; ++iz)
	  {
	    auto z = iz * 0.01;
	    auto h1s = __gnu_cxx::cyl_hankel_1(nu, z);
	    auto h1b = beast::cyl_hankel_1(nu, z);
	    std::cout << ' ' << std::setw(width) << nu
		      << ' ' << std::setw(width) << z
		      << ' ' << std::setw(2 * width) << h1s
		      << ' ' << std::setw(2 * width) << h1b
		      << ' ' << std::setw(width) << std::abs((h1s - h1b) / h1b)
		      << std::endl;
	  }
      }
  }

int
main()
{
  RunCylHankel1<double>();
}
