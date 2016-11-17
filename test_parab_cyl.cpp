/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_parab_cyl test_parab_cyl.cpp
./test_parab_cyl > test_parab_cyl.txt

*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include "sf_parab_cyl.tcc"

template<typename _Tp>
  void
  test_parab_cyl(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__max_digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    // _Tp a = 0.5; // blows up a factor!
    _Tp a = 1.25;
    for (int i = 0; i <= +200; ++i)
      {
	auto z = _Tp{0.1L} * i;
	auto [U, V] = std::__detail::__parabolic_cylinder(a, z);
	std::cout << ' ' << std::setw(width) << z
		  << ' ' << std::setw(width) << U
		  << ' ' << std::setw(width) << V;
	std::cout << '\n';
      }
  }

int
main()
{
  test_parab_cyl(1.0);
}
