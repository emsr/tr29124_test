/**
 *
 */

/*
These are from cephes.
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <test_struve_old.h>

#include <emsr/sf_bessel.h>
#include <emsr/numeric_limits.h>

#include "hyperg_1F2.cpp"
#include "hyperg_2F0.cpp"
#include "hyperg_3F0.cpp"

double PI = 3.141592653589793238462643383279502884195;

double
yv(double nu, double x)
{ return emsr::cyl_neumann(nu, x); }

double
jv(double nu, double x)
{ return emsr::cyl_bessel_j(nu, x); }

double
struve(double nu, double x)
{
  double y, ya;
  double onef2err, threef0err;

  double f = std::floor(nu);
  if (nu < 0.0 && nu - f == 0.5)
    {
      y = jv(-nu, x);
      f = 1.0 - f;
      auto g = 2.0 * std::floor(f / 2.0);
      if (g != f)
        y = -y;
      return y;
    }
  double t = 0.25 * x * x;
  f = std::abs(x);
  double g = 1.5 * std::abs(nu);
  if (f > 30.0 && f > g)
    {
      onef2err = 1.0e38;
      y = 0.0;
    }
  else
    y = hyperg_1f2(1.0, 1.5, 1.5 + nu, -t, onef2err);

  if (f < 18.0 || x < 0.0)
    {
      threef0err = 1.0e38;
      ya = 0.0;
    }
  else
    ya = hyperg_3f0(1.0, 0.5, 0.5 - nu, -1.0 / t, threef0err);

  f = std::sqrt(PI);
  double h = std::pow(0.5 * x, nu - 1.0);

  if (onef2err <= threef0err)
    {
      g = gamma(nu + 1.5);
      y *= h * t / ( 0.5 * f * g );
      return y;
    }
  else
    {
      g = gamma(nu + 0.5);
      ya *= h / (f * g);
      ya += yv(nu, x);
      return ya;
    }
}

int
main()
{
  using Tp = double;

  std::cout.precision(emsr::digits10<Tp>());
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  std::vector<Tp> nuvec
  {
    Tp{0},
    Tp{1} / Tp{3},
    Tp{1} / Tp{2},
    Tp{2} / Tp{3},
    Tp{1},
    Tp{2},
    Tp{5},
    Tp{10},
    Tp{20},
    Tp{50},
    Tp{100},
  };
  for (auto nu : nuvec)
    {
      std::cout << "\n  nu = " << nu << '\n';
      std::cout << ' ' << std::setw(width) << "x";
      std::cout << ' ' << std::setw(width) << "H_nu(x)";
      std::cout << '\n';
      for (unsigned int i = 0; i <= 1000; ++i)
        {
          auto x = Tp(i) / Tp{100};
          std::cout << ' ' << std::setw(width) << x;
          try
            {
              auto h = struve(nu, x);
              std::cout << ' ' << std::setw(width) << h;
            }
          catch (std::exception& err)
            {
              std::cerr << err.what() << '\n';
            }
          std::cout << '\n';
        }
    }
}

