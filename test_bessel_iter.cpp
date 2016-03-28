// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_bessel_iter test_bessel_iter.cpp

// ./test_bessel_iter > test_bessel_iter.txt

// g++ -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_bessel_iter test_bessel_iter.cpp

// ./test_bessel_iter > test_bessel_iter.txt

#include <iostream>
#include <iomanip>
#include "polynomial.h"

// This thing will develop a set of rational functions to express
// Cylinder functions of order nu in terms to rational functions of nu and 2/z
// multiplying Cylinder functions of order nu + n and nu + n + 1.

/*
template<typename _Tp>
  std::pair<, >
  bumpup(int k, )
  {
    __gnu_cxx::_Polynomial<__gnu_cxx::_Polynomial<_Tp>>
    thing
    {
      {{0}, {_Tp(k + 1), _Tp{1}}},
      {{-1}}
    };
  }
*/

int
main()
{
  using _Tp = double;
  using _RatTp = __gnu_cxx::_Polynomial<__gnu_cxx::_Polynomial<_Tp>>;

  std::vector<_RatTp> thing;
  thing.push_back({{0}, {_Tp{1}, _Tp{1}}});
  thing.push_back({{-1}});
  for (int __k = 2; __k <= 25; ++__k)
    {
      _RatTp P0{{0}, {_Tp(__k), _Tp{1}}};
      thing.push_back(P0 * thing[__k - 2] + thing[__k - 1]);
      _RatTp P1{{-1}};
      thing.push_back(P1 * thing[__k - 2]);
    }
  for (int __k = 0; __k < 25; ++__k)
    {
      std::cout << "nu + " << std::setw(2) << __k + 1 << ": " << thing[2 * __k] << '\n';
      std::cout << "nu + " << std::setw(2) << __k + 2 << ": " << thing[2 * __k + 1] << '\n';
      std::cout << '\n';
    }
}
