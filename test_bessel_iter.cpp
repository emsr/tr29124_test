// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_bessel_iter test_bessel_iter.cpp

// ./test_bessel_iter > test_bessel_iter.txt

// g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_bessel_iter test_bessel_iter.cpp

// ./test_bessel_iter > test_bessel_iter.txt

#include <iostream>
#include <iomanip>
#include "polynomial.h"

// This thing1 will develop a set of rational functions to express
// Cylinder functions of order nu in terms to rational functions of nu and 2/z
// multiplying Cylinder functions of order nu + n and nu + n + 1.


int
main()
{
  using _Tp = double;
  using _RatTp = __gnu_cxx::_Polynomial<__gnu_cxx::_Polynomial<_Tp>>;

  std::vector<_RatTp> thing1;
  thing1.push_back({{0}, {_Tp{1}, _Tp{1}}});
  thing1.push_back({{-1}});
  for (int __k = 2; __k <= 25; ++__k)
    {
      _RatTp P0{{0}, {_Tp(__k), _Tp{1}}};
      thing1.push_back(P0 * thing1[__k - 2] + thing1[__k - 1]);
      _RatTp P1{{-1}};
      thing1.push_back(P1 * thing1[__k - 2]);
    }
  for (int __k = 0; __k < 25; ++__k)
    {
      std::cout << "nu + " << std::setw(2) << __k + 1 << ": " << thing1[2 * __k] << '\n';
      std::cout << "nu + " << std::setw(2) << __k + 2 << ": " << thing1[2 * __k + 1] << '\n';
      std::cout << '\n';
    }

  std::vector<_RatTp> thing2;
  thing2.push_back({{-1}});
  thing2.push_back({{0}, {_Tp{1}, _Tp{1}}});
  for (int __k = 2; __k <= 25; ++__k)
    {
      _RatTp Q0{{0}, {_Tp(__k), _Tp{1}}};
      thing2.push_back(Q0 * thing2[__k - 2] - thing2[__k - 1]);
      _RatTp Q1{{-1}};
      thing2.push_back(Q1 * thing2[__k - 2]);
    }
  for (int __k = 0; __k < 25; ++__k)
    {
      std::cout << "nu + " << std::setw(2) << 0 << ": " << thing2[2 * __k] << '\n';
      std::cout << "nu + " << std::setw(2) << 1 << ": " << thing2[2 * __k + 1] << '\n';
      std::cout << '\n';
    }
}
