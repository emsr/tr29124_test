/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hankel_real_arg test_hankel_real_arg.cpp wrap_boost.cpp
./test_hankel_real_arg > test_hankel_real_arg.txt

*/

#include <iostream>
#include <iomanip>
#include <bits/specfun.h>
#include "wrap_boost.h"

int
main()
{
  std::cout << '\n';
  for (auto nu : {-5.0, -2.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0})
    for (auto x : {-5.0, -4.0, -3.0, -2.0, -1.0,/* 0.0,*/ 1.0, 2.0, 3.0, 4.0, 5.0})
      {
	auto H1_boost = beast::cyl_hankel_1(nu, x);
	auto H1_gnu = __gnu_cxx::cyl_hankel_1(nu, x);
	std::cout << ' ' << std::setw(4) << nu
		  << ' ' << std::setw(4) << x
		  << ' ' << std::setw(24) << H1_gnu
		  << ' ' << std::setw(24) << H1_boost
		  << ' ' << std::setw(10) << std::abs(H1_gnu - H1_boost)
		  << '\n';
      }

  std::cout << '\n';
  for (auto n : {0, 1, 2, 5})
    for (auto x : {-5.0, -4.0, -3.0, -2.0, -1.0,/* 0.0,*/ 1.0, 2.0, 3.0, 4.0, 5.0})
      {
	auto h1_boost = beast::sph_hankel_1(n, x);
	auto h1_gnu = __gnu_cxx::sph_hankel_1(n, x);
	std::cout << ' ' << std::setw(4) << n
		  << ' ' << std::setw(4) << x
		  << ' ' << std::setw(24) << h1_gnu
		  << ' ' << std::setw(24) << h1_boost
		  << ' ' << std::setw(10) << std::abs(h1_gnu - h1_boost)
		  << '\n';
      }
}
