/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bessel_iter test_bessel_iter.cpp -lquadmath
./test_bessel_iter > test_bessel_iter.txt

g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_bessel_iter test_bessel_iter.cpp
./test_bessel_iter > test_bessel_iter.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ext/polynomial.h>

// This thing1 will develop a set of rational functions to express
// Cylinder functions of order nu in terms to rational functions of nu and 2/z
// multiplying Cylinder functions of order nu + n and nu + n + 1.


int
main()
{
  using _Tp = double;
  using _RatTp = __gnu_cxx::_Polynomial<__gnu_cxx::_Polynomial<_Tp>>;

  std::cout.precision(std::numeric_limits<_Tp>::digits10);
  auto w = 8 + std::cout.precision();

  std::vector<_RatTp> thing1;
  thing1.push_back({{0}, {_Tp{1}, _Tp{1}}});
  thing1.push_back({{-1}});
  for (int __n = 2; __n <= 25; ++__n)
    {
      _RatTp P0{{0}, {_Tp(__n + 1), _Tp{1}}};
      thing1.push_back(P0 * thing1[__n - 2] + thing1[__n - 1]);
      _RatTp P1{{-1}};
      thing1.push_back(P1 * thing1[__n - 2]);
    }
  for (int __n = 1; __n <= 25; ++__n)
    {
      std::cout << "k = " << 2 * (__n - 1) << ": " << "nu + " << std::setw(2) << __n << ": " << thing1[2 * (__n - 1)] << '\n';
      std::cout << "k = " << 2 * (__n - 1) + 1 << ": " << "nu + " << std::setw(2) << __n + 1 << ": " << thing1[2 * (__n - 1) + 1] << '\n';
      std::cout << '\n';
    }

  for (int __i = 1; __i <= 500; ++__i)
    {
      auto x = 0.1 * __i;
      auto J24 = std::cyl_bessel_j(24.0, x);
      auto J25 = std::cyl_bessel_j(25.0, x);
      auto J0 = std::cyl_bessel_j(0.0, x);
      auto J1 = std::cyl_bessel_j(1.0, x);
      auto J2 = std::cyl_bessel_j(2.0, x);
      auto ladder_J0_J1_J2 = thing1[0](2.0 / x)(0.0) * J1
			   + thing1[1](2.0 / x)(0.0) * J2;
      auto ladder_J0_J24_J25 = thing1[46](2.0 / x)(0.0) * J24
			     + thing1[47](2.0 / x)(0.0) * J25;
      auto ladder_J1_J24_J25 = thing1[46](2.0 / x)(1.0) * J24
			     + thing1[47](2.0 / x)(1.0) * J25;
      std::cout << std::setw(w) << x
		<< std::setw(w) << J0
		<< std::setw(w) << ladder_J0_J1_J2
		<< std::setw(w) << ladder_J0_J24_J25
		<< std::setw(w) << J1
		<< std::setw(w) << ladder_J1_J24_J25
		<< '\n';
    }

  std::vector<_RatTp> thing2;
  thing2.push_back({{-1}});
  thing2.push_back({{0}, {_Tp{1}, _Tp{1}}});
  for (int __n = 2; __n <= 25; ++__n)
    {
      _RatTp Q0{{0}, {_Tp(__n), _Tp{1}}};
      thing2.push_back(Q0 * thing2[__n - 2] - thing2[__n - 1]);
      _RatTp Q1{{-1}};
      thing2.push_back(Q1 * thing2[__n - 2]);
    }
  for (int __n = 0; __n < 25; ++__n)
    {
      std::cout << "C_{nu + " << std::setw(2) << __n + 2 << "} =\n";
      std::cout << "  C_{nu + " << std::setw(2) << 0 << "} * " << thing2[2 * __n] << '\n';
      std::cout << "  C_{nu + " << std::setw(2) << 1 << "} * " << thing2[2 * __n + 1] << '\n';
      std::cout << '\n';
    }
/*
  for (int __i = 1; __i <= 100; ++__i)
    {
      auto x = 0.1 * __i;
      auto J25 = std::cyl_bessel_j(25.0, x);
      auto J0 = std::cyl_bessel_j(0.0, x);
      auto J1 = std::cyl_bessel_j(1.0, x);
    }
*/
}
