#include <numbers>
#include <iostream>

#include <emsr/special_functions.h>

int
main()
{
  const double pi_2 = std::numbers::pi / 2;
  int num_points = 10;
  double a = -1.0;
  double b = +1.0;
  double x;
  for (int i = 0; i < num_points; ++i)
  {
    if (i < num_points / 2)
    {
      auto hth = pi_2 * i / (num_points - 1);
      auto sinhth = std::sin(hth);
      x = a + (b - a) * sinhth * sinhth;
    }
    else
    {
      auto hth = pi_2 * (num_points - i - 1) / (num_points - 1);
      auto sinhth = std::sin(hth);
      x = b - (b - a) * sinhth * sinhth;
    }
    std::cout << '\t'  << x << '\t' << emsr::chebyshev_t(num_points - 1, x) << '\n';
  }
}
