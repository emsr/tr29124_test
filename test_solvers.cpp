/*
$HOME/bin/bin/g++ -std=gnu++17 -g -I. -o test_solvers test_solvers.cpp -lquadmath
./test_solvers > test_solvers.txt
*/

#include "solver.h"
#include <iostream>
#include <iomanip>
#include <bits/numeric_limits.h>

template<typename _Real, unsigned int _Dim>
  bool
  try_solution(std::array<_Real, _Dim> coef, solution_t<_Real> x)
  {
    const auto _S_eps = __gnu_cxx::__epsilon<_Real>();
    if (x.index() == 0)
      return true;
    else if (x.index() == 1)
      {
	_Real y = coef[_Dim - 1];
	for (unsigned int i = _Dim - 1; i > 0; --i)
	  y += coef[i - 1] + x * y;
	return std::abs(y) < _S_eps;
      }
    else if (x.index() == 2)
      {
	std::complex<_Real> y = coef[_Dim - 1];
      }
  }

int
main()
{
  std::cout << '\n';
  auto quad = quadratic(std::experimental::make_array(-2.0, 1.0, 1.0));
  for (int i = 0; i < 2; ++i)
    {
      std::cout << quad[i] << '\n';
    }

  std::cout << '\n';
  auto linquad = quadratic(std::experimental::make_array(-2.0, 1.0, 0.0));
  for (int i = 0; i < 2; ++i)
    {
      std::cout << linquad[i] << '\n';
    }

  std::cout << '\n';
  auto cub = cubic(std::experimental::make_array(24.0, -22.0, -4.0, 2.0));
  for (int i = 0; i < 3; ++i)
    {
      std::cout << cub[i] << '\n';
    }

  std::cout << '\n';
  auto quadcub = cubic(std::experimental::make_array(24.0, -22.0, -4.0, 0.0));
  for (int i = 0; i < 3; ++i)
    {
      std::cout << quadcub[i] << '\n';
    }

  std::cout << '\n';
  auto lincub = cubic(std::experimental::make_array(24.0, -22.0, 0.0, 0.0));
  for (int i = 0; i < 3; ++i)
    {
      std::cout << lincub[i] << '\n';
    }

  std::cout << '\n';
  auto cub2 = cubic(std::experimental::make_array(-4.0, -15.0, 0.0, 1.0));
  for (int i = 0; i < 3; ++i)
    {
      std::cout << cub2[i] << '\n';
    }
}
