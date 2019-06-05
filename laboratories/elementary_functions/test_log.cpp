/**
 *
 */
// Got these from Burkhardt.
// They are used in cdflib.
// These need to be made type generic.
// All these things like expm1, log1p essentially pick off leading coefficients.

#include <cmath>
#include <iostream>
#include <iomanip>

/**
 * Return @f$ log(1 + x) - x @f$.
 * Translate Burkhardt -rlog1.
 */
double
log1pm(double x)
{
  constexpr double a = 0.566749439387324e-01;
  constexpr double b = 0.456512608815524e-01;
  constexpr double p0 = 0.333333333333333e+00;
  constexpr double p1 = -0.224696413112536e+00;
  constexpr double p2 = 0.620886815375787e-02;
  constexpr double q1 = -0.127408923933623e+01;
  constexpr double q2 = 0.354508718369557e+00;

  if (x < -0.39 || x > 0.57)
    return std::log1p(x) - x;
  else
    {
      double h, w1;
      if (x < -0.18)
	{
	  h = (x + 0.3) / 0.7;
	  w1 = a - h * 0.3;
	}
      else if (x > 0.18)
	{
	  h = 0.75 * x - 0.25;
	  w1 = b + h / 3.0;
	}
      else
	{
	  h = x;
	  w1 = 0.0;
	}
      const auto r = h / (h + 2.0);
      const auto t = r * r;
      const auto w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0);
      return -(2.0 * t * (1.0 / (1.0 - r) - r * w) + w1);
    }
}

/**
 * Return @f$ log(x) + 1 - x @f$.
 * Translate Burkhardt -rlog.
 */
double
logp1m(double x)
{
  constexpr double a = 0.566749439387324e-01;
  constexpr double b = 0.456512608815524e-01;
  constexpr double p0 = 0.333333333333333e+00;
  constexpr double p1 = -0.224696413112536e+00;
  constexpr double p2 = 0.620886815375787e-02;
  constexpr double q1 = -0.127408923933623e+01;
  constexpr double q2 = 0.354508718369557e+00;

  if (x < 0.61 || x > 1.57)
    return std::log(x) + 1.0 - x;
  else
    {
      double u, w1;
      if (x < 0.82)
	{
	  u = (x - 0.7) / 0.7;
	  w1 = a - u * 0.3;
	}
      else if (x > 1.18)
	{
	  u = 0.75 * x - 1.0;
	  w1 = b + u / 3.0;
	}
      else
	{
	  u = x - 1.0;
	  w1 = 0.0;
	}
      const auto r = u / (u + 2.0);
      const auto t = r * r;
      const auto w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0);
      return -(2.0 * t * (1.0 / (1.0 - r) - r * w) + w1);
    }
}

int
main()
{
  std::cout.precision(16);
  int w = 8 + std::cout.precision();

  std::cout << "\n\n";
  for (int i = -100; i <= +100; ++i)
    {
      auto x = double(0.01L * i);
      std::cout << ' ' << std::setw(w) << x
		<< ' ' << std::setw(w) << std::log(1.0 + x) - x
		<< ' ' << std::setw(w) << std::log1p(x) - x
		<< ' ' << std::setw(w) << log1pm(x)
		<< '\n';
    }

  std::cout << "\n\n";
  for (int i = -100; i <= +100; ++i)
    {
      auto x = double(1.0L + 0.01L * i);
      std::cout << ' ' << std::setw(w) << x
		<< ' ' << std::setw(w) << std::log(x) + 1 - x
		<< ' ' << std::setw(w) << logp1m(x)
		<< '\n';
    }
}
