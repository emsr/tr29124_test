/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I.. -o test_fft test_fft.cpp -lquadmath
./test_fft > test_fft.txt
*/

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include "fft.h"

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10<double>());
  auto w = 8 + std::cout.precision();
  auto cw = 4 + 2 * w;

  std::default_random_engine re;
  std::uniform_real_distribution<double> ud(0, 2 * M_PI);
  auto gen = [&ud, &re]()->double{ return ud(re); };
  std::vector<std::complex<double>> vec;
  vec.reserve(1000);
  for (int i = 0; i < 1000; ++i)
    vec.push_back(std::polar(1.0, gen()));

  auto xform = vec;
  __gnu_cxx::fft(xform);

  auto iform = xform;
  __gnu_cxx::ifft(iform);

  for (int i = 0; i < 1000; ++i)
    std::cout << ' ' << std::setw(cw) << vec[i]
	      << ' ' << std::setw(cw) << xform[i]
	      << ' ' << std::setw(cw) << iform[i]
	      << ' ' << std::setw(cw) << vec[i] - iform[i]
	      << ' ' << std::setw(w) << std::abs(vec[i] - iform[i])
	      << '\n';
}
