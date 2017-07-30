/*
$HOME/bin/bin/g++ -std=gnu++17 -g -I. -o test_solvers test_solvers.cpp -lquadmath
./test_solvers > test_solvers.txt
*/

#include "solver.h"
#include <iostream>
#include <iomanip>
#include <bits/numeric_limits.h>

template<typename _Real, unsigned long int _Dim>
  bool
  try_solution(std::array<_Real, _Dim> coef, const solution_t<_Real>& z)
  {
    const auto _S_eps = __gnu_cxx::__epsilon<_Real>();
    if (z.index() == 0)
      return true;
    else if (z.index() == 1)
      {
	const _Real x = std::get<1>(z);
	_Real y = coef[_Dim - 1];
	for (unsigned int i = _Dim - 1; i > 0; --i)
	  y = coef[i - 1] + x * y;
//std::cout << ' ' << y << '\n';
	return std::abs(y) < _S_eps;
      }
    else if (z.index() == 2)
      {
	const std::complex<_Real> x = std::get<2>(z);
	std::complex<_Real> y = coef[_Dim - 1];
	for (unsigned int i = _Dim - 1; i > 0; --i)
	  y = coef[i - 1] + x * y;
	return std::abs(y) < _S_eps;
//std::cout << ' ' << y << '\n';
      }
  }

int
main()
{
  using std::experimental::make_array;

  std::cout << std::boolalpha;
  std::cout.precision(__gnu_cxx::__digits10<double>());
  std::cout << std::showpoint << std::scientific;
  auto w = 8 + std::cout.precision();

  std::cout << '\n';
  const auto quad_coef = make_array(-2.0, 1.0, 1.0);
  auto quad = quadratic(quad_coef);
  for (int i = 0; i < 2; ++i)
    {
      bool ok = try_solution(quad_coef, quad[i]);
      std::cout << std::setw(w) << quad[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto cquad_coef = make_array(2.0, 1.0, 1.0);
  auto cquad = quadratic(cquad_coef);
  for (int i = 0; i < 2; ++i)
    {
      bool ok = try_solution(cquad_coef, cquad[i]);
      std::cout << std::setw(w) << cquad[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto linquad_coef = make_array(-2.0, 1.0, 0.0);
  auto linquad = quadratic(linquad_coef);
  for (int i = 0; i < 2; ++i)
    {
      bool ok = try_solution(linquad_coef, linquad[i]);
      std::cout << std::setw(w) << linquad[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto cub_coef = make_array(24.0, -22.0, -4.0, 2.0);
  auto cub = cubic(cub_coef);
  for (int i = 0; i < 3; ++i)
    {
      bool ok = try_solution(cub_coef, cub[i]);
      std::cout << std::setw(w) << cub[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto quadcub_coef = make_array(24.0, -22.0, -4.0, 0.0);
  auto quadcub = cubic(quadcub_coef);
  for (int i = 0; i < 3; ++i)
    {
      bool ok = try_solution(quadcub_coef, quadcub[i]);
      std::cout << std::setw(w) << quadcub[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto lincub_coef = make_array(24.0, -22.0, 0.0, 0.0);
  auto lincub = cubic(lincub_coef);
  for (int i = 0; i < 3; ++i)
    {
      bool ok = try_solution(lincub_coef, lincub[i]);
      std::cout << std::setw(w) << lincub[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto cub2_coef = make_array(-4.0, -15.0, 0.0, 1.0);
  auto cub2 = cubic(cub2_coef);
  for (int i = 0; i < 3; ++i)
    {
      bool ok = try_solution(cub2_coef, cub2[i]);
      std::cout << std::setw(w) << cub2[i] << "  good: " << ok << '\n';
    }
/*
  std::cout << '\n';
  const auto quart_coef = make_array(-30.0, 31.0, 5.0, -7.0, 1.0);
  auto quart = quartic(quart_coef);
  for (int i = 0; i < 4; ++i)
    {
      bool ok = try_solution(quart_coef, quart[i]);
      std::cout << std::setw(w) << quart[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto cubquart_coef = make_array(-30.0, 31.0, 5.0, -7.0, 0.0);
  auto cubquart = quartic(cubquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      bool ok = try_solution(cubquart_coef, cubquart[i]);
      std::cout << std::setw(w) << cubquart[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto quadquart_coef = make_array(-30.0, 31.0, 5.0, 0.0, 0.0);
  auto quadquart = quartic(quadquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      bool ok = try_solution(quadquart_coef, quadquart[i]);
      std::cout << std::setw(w) << quadquart[i] << "  good: " << ok << '\n';
    }

  std::cout << '\n';
  const auto linquart_coef = make_array(-30.0, 31.0, 0.0, 0.0, 0.0;
  auto linquart = quartic(linquart_coef));
  for (int i = 0; i < 4; ++i)
    {
      bool ok = try_solution(linquart_coef, linquart[i]);
      std::cout << std::setw(w) << linquart[i] << "  good: " << ok << '\n';
    }
*/
}
