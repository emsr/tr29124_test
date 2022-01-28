/**
 *
 */

#include <limits>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/sf_gamma.h>

#include <wrap_boost.h>

namespace emsr
{

  template<typename Tp>
    Tp
    rising_factorial_upper_prod(int a, int n)
    {
      if (a < n)
	return Tp{0};
      else
	{
	  auto prod = 1;
	  for (int k = 0; k < n; ++k)
	    prod *= a--;
	  return prod;
	}
    }

  template<typename Tp>
    Tp
    rising_factorial_prod(Tp a, int n)
    {
      auto prod = Tp{1};
      for (int k = 0; k < n; ++k)
	prod *= a++;
      return prod;
    }

  template<typename Tp>
    Tp
    rising_factorial_fake(Tp a, Tp x)
    {
      auto n = int(std::nearbyint(x));
      if (Tp(n) == x)
	{
	  if (n == 0)
	    return Tp{1};
	  else
	    {
	      auto m = int(std::nearbyint(a));
	      if (int(m) == a)
		return rising_factorial_prod<Tp>(m, n);
	      else
		return rising_factorial_prod(a, n);
	    }
	}
      else
	return emsr::detail::gamma(a + x)
	     / emsr::detail::gamma(a);
    }

}

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 8 + std::cout.precision();

  std::vector<double> xv{0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};

  // Test call (log form).
  std::cout << "\nLogarithmic form\n";
  for (int ia = 0; ia <= +500; ++ia)
    {
      auto a = ia * 0.01;
      std::cout << '\n';
      for (auto x : xv)
	{
	  auto pochg = emsr::rising_factorial(a, x);
	  auto pochb = beast::rising_factorial(a, x);
	  std::cout << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << x
		    << ' ' << std::setw(width) << pochg
		    << ' ' << std::setw(width) << pochb
		    << ' ' << std::setw(width) << pochg - pochb
		    << ' ' << std::setw(width) << (pochg - pochb) / pochb
		    << '\n';
	}
    }

  // Test naive product form.
  std::cout << "\nProduct form\n";
  for (int ia = 0; ia <= +500; ++ia)
    {
      auto a = ia * 0.01;
      std::cout << '\n';
      for (auto x : xv)
	{
	  auto pochg = emsr::rising_factorial_fake(a, x);
	  auto pochb = beast::rising_factorial(a, x);
	  std::cout << ' ' << std::setw(width) << a
		    << ' ' << std::setw(width) << x
		    << ' ' << std::setw(width) << pochg
		    << ' ' << std::setw(width) << pochb
		    << ' ' << std::setw(width) << pochg - pochb
		    << ' ' << std::setw(width) << (pochg - pochb) / pochb
		    << '\n';
	}
    }
}
