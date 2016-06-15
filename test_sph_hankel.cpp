/*
$HOME/bin_tr29124/bin/g++ -g -o test_sph_hankel test_sph_hankel.cpp -lquadmath boost_wrap.cpp
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_sph_hankel > test_sph_hankel.txt

g++ -std=c++14 -g -o test_sph_hankel test_sph_hankel.cpp -lquadmath
./test_sph_hankel > test_sph_hankel.txt
*/

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include "boost_wrap.h"

template<typename _Tp>
  void
  RunSphHankel1()
  {
    std::cout.precision(std::numeric_limits<double>::digits10);
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    auto max_max_abs_frac = _Tp{-1};
    for (int inu = 0; inu <= +500; ++inu)
      {
	auto max_abs_frac = _Tp{-1};
	auto nu = inu * 0.1;
	std::cout << '\n';
	for (int iz = 1; iz <= +1000; ++iz)
	  {
	    auto z = iz * 0.01;
	    auto h1s = __gnu_cxx::sph_hankel_1(nu, z);
	    auto h1b = beast::sph_hankel_1(nu, z);
	    auto abs_frac = std::abs((h1s - h1b) / h1b);
	    std::cout << ' ' << std::setw(width) << nu
		      << ' ' << std::setw(width) << z
		      << ' ' << std::setw(2 * width) << h1s
		      << ' ' << std::setw(2 * width) << h1b
		      << ' ' << std::setw(width) << abs_frac
		      << std::endl;
	    if (abs_frac > max_abs_frac)
	      max_abs_frac = abs_frac;
	  }
	std::cout << " max_abs_frac = " << std::setw(width) << max_abs_frac << '\n';
	if (max_abs_frac > max_max_abs_frac)
	  max_max_abs_frac = max_abs_frac;
      }
    std::cout << " max_max_abs_frac = " << std::setw(width) << max_max_abs_frac << '\n';
  }

int
main()
{
  RunSphHankel1<double>();
}
